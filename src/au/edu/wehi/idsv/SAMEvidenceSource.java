package au.edu.wehi.idsv;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.List;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;

import com.google.common.collect.Lists;

/**
 * Structural variation evidence based on read pairs from a single SAM/BAM.  
 * @author cameron.d
 *
 */
public class SAMEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(SAMEvidenceSource.class);
	public static final String FLAG_NAME = "extracted";
	private final ProcessingContext processContext;
	private final File input;
	private final boolean isTumour;
	private RelevantMetrics metrics;
	public SAMEvidenceSource(ProcessingContext processContext, File file, boolean isTumour) {
		super(processContext, file);
		this.processContext = processContext;
		this.input = file;
		this.isTumour = isTumour;
	}
	/**
	 * Ensures that all structural variation evidence has been extracted from the input file 
	 */
	public void ensureEvidenceExtracted() {
		if (!isProcessingComplete()) {
			process();
		} else {
			log.debug("No extraction required for ", input);
		}
	}
	@Override
	protected void process() {
		try (ExtractEvidence extract = new ExtractEvidence()) {
			log.info("START extract evidence for ", input);
			extract.extractEvidence();
			log.info("SUCCESS extract evidence for ", input);
		}
	}
	@Override
	protected boolean isProcessingComplete() {
		boolean done = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		done &= checkIntermediate(fsc.getMetrics(input), input);
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				done &= checkIntermediate(fsc.getSVBamForChr(input, seq.getSequenceName()), input);
				done &= checkIntermediate(fsc.getMateBamForChr(input, seq.getSequenceName()), input);
				done &= checkIntermediate(fsc.getRealignmentFastqForChr(input, seq.getSequenceName()), input);
			}
		} else {
			done &= checkIntermediate(fsc.getSVBam(input), input);
			done &= checkIntermediate(fsc.getMateBam(input), input);
			done &= checkIntermediate(fsc.getRealignmentFastq(input), input);
		}
		done &= checkIntermediate(fsc.getFlagFile(input, FLAG_NAME), input);
		return done;
	}
	public RelevantMetrics getMetrics() {
		if (metrics == null) {
			ensureEvidenceExtracted();
			metrics = new PicardMetrics(processContext.getFileSystemContext().getMetrics(input));
		}
		return metrics;
	}
	public File getSourceFile() { return input; }
	public boolean isTumour() { return isTumour; }
	@Override
	protected DirectedEvidenceFileIterator perChrIterator(String chr) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return new DirectedEvidenceFileIterator(processContext, this,
				fsc.getSVBamForChr(input, chr),
				fsc.getMateBamForChr(input, chr),
				isRealignmentComplete() ? fsc.getRealignmentBamForChr(input, chr) : null,
				null);
	}
	@Override
	protected DirectedEvidenceFileIterator singleFileIterator() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return new DirectedEvidenceFileIterator(processContext, this,
				fsc.getSVBam(input),
				fsc.getMateBam(input),
				isRealignmentComplete() ? fsc.getRealignmentBam(input) : null,
				null);
	}
	/**
	 * Extracts reads supporting structural variation into intermediate files.
	 * By sorting DP and OEA read pairs by the coordinate of the mapped mate,
	 * putative directed breakpoints can be assembled by downstream processes
	 * in a single pass of the intermediate files.
	 * @author Daniel Cameron
	 *
	 */
	private class ExtractEvidence implements Closeable {
		private SamReader reader = null;
		private ReferenceSequenceFileWalker referenceWalker = null;
		private final List<SAMFileWriter> writers = Lists.newArrayList();
		private final List<SAMFileWriter> matewriters = Lists.newArrayList();
		private final List<SortingCollection<SAMRecord>> mateRecordBuffer = Lists.newArrayList();
		private final List<FastqBreakpointWriter> realignmentWriters = Lists.newArrayList();
		public void close() {
			tryClose(reader);
			reader = null;
			tryClose(referenceWalker);
			referenceWalker = null;
			for (SAMFileWriter w : writers) {
				tryClose(w);
			}
			writers.clear();
	    	for (SAMFileWriter w : matewriters) {
				tryClose(w);
			}
	    	matewriters.clear();
	    	mateRecordBuffer.clear();
	    	for (FastqBreakpointWriter w : realignmentWriters) {
				tryClose(w);
			}
	    	realignmentWriters.clear();
		}
		private void tryClose(Closeable toClose) {
			try {
				if (toClose != null) toClose.close();
	    	} catch (IOException e) {
	    		log.equals(e);
	    	}
		}
		/**
		 * Deletes all output files
		 * Deleting output files if there is an error prevents downstream issues
		 * with partially written intermediate files 
		 */
		private void deleteOutput() {
			close(); // close any file handles that are still around
			FileSystemContext fsc = processContext.getFileSystemContext();
			tryDelete(fsc.getFlagFile(input, FLAG_NAME));
			tryDelete(fsc.getMetrics(input));
			tryDelete(fsc.getSVBam(input));
			tryDelete(fsc.getMateBam(input));
			tryDelete(fsc.getRealignmentFastq(input));
			for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
				String seq = seqr.getSequenceName();
				tryDelete(fsc.getSVBamForChr(input, seq));
				tryDelete(fsc.getMateBamForChr(input, seq));
				tryDelete(fsc.getRealignmentFastqForChr(input, seq));
			}
		}
		private void tryDelete(File f) {
			try {
				if (f.exists()) {
					if (!f.delete()) {
						log.error("Unable to delete intermediate file ", f,  " during rollback.");
					}
				}
			} catch (Exception e) {
				log.error(e, "Unable to delete intermediate file ", f,  " during rollback.");
			}
		}
	    private void writeToFile(SortingCollection<SAMRecord> sc, SAMFileWriter writer) {
	    	sc.doneAdding();
	    	final CloseableIterator<SAMRecord> it = sc.iterator();
			while (it.hasNext()) {
				writer.addAlignment(it.next());
			}
	    }
		public void extractEvidence() {
	    	try {
	    		deleteOutput(); // trash any left-over files
	    		
		    	IOUtil.assertFileIsWritable(processContext.getFileSystemContext().getMetrics(input));
		    	
		    	reader = processContext.getSamReaderFactory().open(input);
		    	final SAMFileHeader header = reader.getFileHeader();
		    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
		    	referenceWalker = new ReferenceSequenceFileWalker(processContext.getReferenceFile());
		    	final InsertSizeMetricsCollector metrics = PicardMetrics.createCollector(header);
		    	
		    	dictionary.assertSameDictionary(processContext.getReference().getSequenceDictionary());
		    	
		    	createOutputWriters(header, dictionary);
		    	
		    	processInputRecords(dictionary, metrics);
		    	
		    	flush();
		    	
		    	writeMetrics(metrics);
		    	
		    	if (!processContext.getFileSystemContext().getFlagFile(input, FLAG_NAME).createNewFile()) {
		    		log.error("Failed to create flag file ", processContext.getFileSystemContext().getFlagFile(input, FLAG_NAME));
		    	}
	    	} catch (Exception e) {
	    		log.error(e, "Unable to extract evidence for ", input);
	    		deleteOutput();
	    		throw new RuntimeException(e);
	    	} finally {
	    		close();
	    	}
	    }
		private void writeMetrics(InsertSizeMetricsCollector metrics) throws IOException {
			log.info("Writing metrics");
			metrics.finish();
			PicardMetrics.save(metrics, processContext.<InsertSizeMetrics, Integer>createMetricsFile(), processContext.getFileSystemContext().getMetrics(input));
		}
		private void flush() throws IOException {
			reader.close();
			reader = null;
			for (SAMFileWriter w : writers) {
				w.close();
			}
			writers.clear();
			for (FastqBreakpointWriter w : realignmentWriters) {
				w.close();
			}
			realignmentWriters.clear();
			log.info("Sorting sv mates");
			for (int i = 0; i < mateRecordBuffer.size(); i++) {
				// TODO: multi-thread this?
				SortingCollection<SAMRecord> buffer = mateRecordBuffer.get(i);
				SAMFileWriter w = matewriters.get(i);
				w.setProgressLogger(new ProgressLogger(log));
				writeToFile(buffer, w);
				w.close();
				matewriters.set(i, null);
				mateRecordBuffer.set(i,  null);
			}
			matewriters.clear();
			mateRecordBuffer.clear();
		}
		private void processInputRecords(final SAMSequenceDictionary dictionary, final InsertSizeMetricsCollector metrics) {
			// output all to a single file, or one per chr + one for unmapped 
			assert(matewriters.size() == mateRecordBuffer.size());
			assert(writers.size() == dictionary.getSequences().size() + 1 || writers.size() == 1);
			assert(matewriters.size() == dictionary.getSequences().size() || matewriters.size() == 1);
			assert(realignmentWriters.size() == dictionary.getSequences().size() || realignmentWriters.size() == 1);
			
			// Traverse the input file
			final ProgressLogger progress = new ProgressLogger(log);
			final SAMRecordIterator iter = reader.iterator().assertSorted(SortOrder.coordinate);
			while (iter.hasNext()) {
				SAMRecord record = iter.next();
				int offset = record.getReadUnmappedFlag() ? dictionary.size() : record.getReferenceIndex();
				boolean sc = SAMRecordUtil.isSoftClipLengthAtLeast(record, processContext.getSoftClipParameters().minLength);				
				boolean badpair = SAMRecordUtil.isPartOfNonReferenceReadPair(record);
				if (sc || badpair) {
					SAMRecordUtil.ensureNmTag(referenceWalker, record);
					writers.get(offset % writers.size()).addAlignment(record);
				}
				if (badpair) {
					if (!record.getMateUnmappedFlag()) {
						mateRecordBuffer.get(record.getMateReferenceIndex() % mateRecordBuffer.size()).add(record);
					}
				}
				if (SAMRecordUtil.getStartSoftClipLength(record) > 0) {
					SoftClipEvidence startEvidence = new SoftClipEvidence(processContext, SAMEvidenceSource.this, BreakendDirection.Backward, record);
					if (processContext.getSoftClipParameters().meetsCritera(startEvidence) && processContext.getRealignmentParameters().shouldRealignBreakend(startEvidence)) {
						realignmentWriters.get(offset % realignmentWriters.size()).write(startEvidence);
					}
				}
				if (SAMRecordUtil.getEndSoftClipLength(record) > 0) {
					SoftClipEvidence endEvidence = new SoftClipEvidence(processContext, SAMEvidenceSource.this, BreakendDirection.Forward, record);
					if (processContext.getSoftClipParameters().meetsCritera(endEvidence) && processContext.getRealignmentParameters().shouldRealignBreakend(endEvidence)) {
						realignmentWriters.get(offset % realignmentWriters.size()).write(endEvidence);
					}
				}
				metrics.acceptRecord(record, null);
				progress.record(record);
			}
		}
		private void createOutputWriters(final SAMFileHeader header, final SAMSequenceDictionary dictionary) throws IOException {
			final SAMFileHeader svHeader = header.clone();
			svHeader.setSortOrder(SortOrder.coordinate);
			final SAMFileHeader mateHeader = header.clone();
			mateHeader.setSortOrder(SortOrder.unsorted);
			if (processContext.shouldProcessPerChromosome()) {
				createOutputWritersPerChromosome(dictionary, svHeader, mateHeader);
			} else {
				createOutputWriterPerGenome(svHeader, mateHeader);
			}
		}
		private void createOutputWriterPerGenome(final SAMFileHeader svHeader, final SAMFileHeader mateHeader) throws IOException {
			// all writers map to the same one
			writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSVBam(input)));
			matewriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(mateHeader, false, processContext.getFileSystemContext().getMateBam(input)));
			mateRecordBuffer.add(SortingCollection.newInstance(
					SAMRecord.class,
					new BAMRecordCodec(mateHeader),
					new SAMRecordMateCoordinateComparator(),
					processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
					processContext.getFileSystemContext().getTemporaryDirectory()));
			realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastq(input))));
		}
		private void createOutputWritersPerChromosome(
				final SAMSequenceDictionary dictionary,
				final SAMFileHeader svHeader, final SAMFileHeader mateHeader)
				throws IOException {
			for (SAMSequenceRecord seq : dictionary.getSequences()) {
				writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSVBamForChr(input, seq.getSequenceName())));
				matewriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(mateHeader, false, processContext.getFileSystemContext().getMateBamForChr(input, seq.getSequenceName())));
				mateRecordBuffer.add(SortingCollection.newInstance(
						SAMRecord.class,
						new BAMRecordCodec(mateHeader),
						new SAMRecordMateCoordinateComparator(),
						processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
						processContext.getFileSystemContext().getTemporaryDirectory()));
				realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastqForChr(input, seq.getSequenceName()))));
			}
			writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSVBamForChr(input, "unmapped")));
		}
	}
}