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
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class SAMRecordEvidenceFile implements Iterable<DirectedEvidence> {
	private static final Log log = Log.getInstance(SAMRecordEvidenceFile.class);
	private final ProcessingContext processContext;
	private final File input;
	private RelevantMetrics metrics;
	public SAMRecordEvidenceFile(ProcessingContext processContext, File file) {
		this.processContext = processContext;
		this.input = file;
	}
	/**
	 * Ensures that all structural variation evidence has been extracted from the input file 
	 */
	public void ensureEvidenceExtracted() {
		if (!extractionHappy()) {
			try (ExtractEvidence extract = new ExtractEvidence()) {
				log.info("Extracting evidence for %s", input);
				extract.extractEvidence();
			}
		}
	}
	private boolean checkHappy(File file, File source) {
		if (!file.exists()) {
			log.debug("Missing intermediate ", file);
			return false;
		}
		if (file.lastModified() < source.lastModified()) {
			log.info(source, " has a more recent timestamp than ", file, ". Considering ", file, " out of date.");
			return false;
		}
		return true;
	}
	private boolean extractionHappy() {
		boolean happy = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		happy |= checkHappy(fsc.getMetrics(input), input);
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				happy |= checkHappy(fsc.getSVBamForChr(input, seq.getSequenceName()), input);
				happy |= checkHappy(fsc.getMateBamForChr(input, seq.getSequenceName()), input);
			}
		} else {
			happy |= checkHappy(fsc.getSVBam(input), input);
			happy |= checkHappy(fsc.getMateBam(input), input);
		}
		return happy;
	}
	private boolean realignmentHappy() {
		boolean happy = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				happy |= checkHappy(fsc.getRealignmentBamForChr(input, seq.getSequenceName()), fsc.getSVBamForChr(input, seq.getSequenceName()));
			}
		} else {
			happy |= checkHappy(fsc.getRealignmentBam(input), fsc.getSVBam(input));
		}
		return happy;
	}
	public RelevantMetrics getMetrics() {
		if (metrics == null) {
			ensureEvidenceExtracted();
			metrics = new RelevantMetrics(processContext.getFileSystemContext().getMetrics(input));
		}
		return metrics;
	}
	public File getSourceFile() { return input; }
	@Override
	public Iterator<DirectedEvidence> iterator() {
		if (!extractionHappy()) {
			throw new IllegalStateException(String.format("Missing intermediate files for ", input));
		} else if (!realignmentHappy()) {
			throw new IllegalStateException(String.format("Missing expected realignment. Have the fastq files in %s been aligned?", processContext.getFileSystemContext().getIntermediateDirectory(input)));
		}
		// TODO Auto-generated method stub
		return null;
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
		private List<SAMFileWriter> writers = Lists.newArrayList();
		private List<SAMFileWriter> matewriters = Lists.newArrayList();
		private List<SortingCollection<SAMRecord>> mateRecordBuffer = Lists.newArrayList();
		/**
		 * Files to delete if there is an error during processing.
		 * Deleting output files if there is an error prevents downstream issues
		 * with partially written intermediate files 
		 */
		private Set<File> deleteOnError = Sets.newHashSet();
		public void close() {
			tryClose(reader);
			reader = null;
			tryClose(referenceWalker);
			referenceWalker = null;
			for (SAMFileWriter w : writers) {
				tryClose(w);
			}
	    	writers = null;
	    	for (SAMFileWriter w : matewriters) {
				tryClose(w);
			}
	    	matewriters = null;
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
		 */
		private void deleteOutput() {
			close(); // close any file handles that are still around
			for (File f : deleteOnError) {
				f = new File(f.getAbsolutePath());
				if (f.exists()) {
					try {
						if (!f.delete()) {
							log.error("Unable to delete intermediate file ", f,  " during rollback.");
						}
					} catch (Exception e) {
						log.error(e, "Unable to delete intermediate file ", f,  " during rollback.");
					}
				}
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
		    	IOUtil.assertFileIsWritable(processContext.getFileSystemContext().getMetrics(input));
		    	
		    	reader = processContext.getSamReaderFactory().open(input);
		    	final SAMFileHeader header = reader.getFileHeader();
		    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
		    	referenceWalker = new ReferenceSequenceFileWalker(processContext.getReferenceFile());
		    	final InsertSizeMetricsCollector metrics = RelevantMetrics.createCollector(header);
		    	
		    	dictionary.assertSameDictionary(processContext.getReference().getSequenceDictionary());
		    	
		    	createSamOutputWriters(header, dictionary);
		    	
		    	processInputRecords(dictionary, metrics);
		    	
		    	flush();
		    	
		    	writeMetrics(metrics);
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
			RelevantMetrics.save(metrics, processContext.<InsertSizeMetrics, Integer>createMetricsFile(), processContext.getFileSystemContext().getMetrics(input));
		}
		private void flush() throws IOException {
			reader.close();
			reader = null;
			for (SAMFileWriter w : writers) {
				w.close();
			}
			writers.clear();
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
			mateRecordBuffer.clear();
			matewriters.clear();
		}
		private void processInputRecords(final SAMSequenceDictionary dictionary, final InsertSizeMetricsCollector metrics) {
			// output all to a single file, or one per chr + one for unmapped 
			assert(matewriters.size() == mateRecordBuffer.size());
			assert(writers.size() == dictionary.getSequences().size() + 1 || writers.size() == 1);
			assert(matewriters.size() == dictionary.getSequences().size() || matewriters.size() == 1);
			
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
				metrics.acceptRecord(record, null);
				progress.record(record);
			}
		}
		private void createSamOutputWriters(final SAMFileHeader header, final SAMSequenceDictionary dictionary) throws IOException {
			final SAMFileHeader svHeader = header.clone();
			svHeader.setSortOrder(SortOrder.coordinate);
			final SAMFileHeader mateHeader = header.clone();
			mateHeader.setSortOrder(SortOrder.unsorted);
			if (processContext.shouldProcessPerChromosome()) {
				createSamOutputWritersPerChromosome(dictionary, svHeader, mateHeader);
			} else {
				createSamOutputWriterPerGenome(svHeader, mateHeader);
			}
		}
		private void createSamOutputWriterPerGenome(final SAMFileHeader svHeader, final SAMFileHeader mateHeader) throws IOException {
			// all writers map to the same one
			writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, track(processContext.getFileSystemContext().getSVBam(input))));
			matewriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(mateHeader, false, track(processContext.getFileSystemContext().getMateBam(input))));
			mateRecordBuffer.add(SortingCollection.newInstance(
					SAMRecord.class,
					new BAMRecordCodec(mateHeader),
					new SAMRecordMateCoordinateComparator(),
					processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
					processContext.getFileSystemContext().getTemporaryDirectory()));
		}
		private void createSamOutputWritersPerChromosome(
				final SAMSequenceDictionary dictionary,
				final SAMFileHeader svHeader, final SAMFileHeader mateHeader)
				throws IOException {
			for (SAMSequenceRecord seq : dictionary.getSequences()) {
				writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, track(processContext.getFileSystemContext().getSVBamForChr(input, seq.getSequenceName()))));
				matewriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(mateHeader, false, track(processContext.getFileSystemContext().getMateBamForChr(input, seq.getSequenceName()))));
				mateRecordBuffer.add(SortingCollection.newInstance(
						SAMRecord.class,
						new BAMRecordCodec(mateHeader),
						new SAMRecordMateCoordinateComparator(),
						processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
						processContext.getFileSystemContext().getTemporaryDirectory()));
			}
			writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, track(processContext.getFileSystemContext().getSVBamForChr(input, "unmapped"))));
		}
		/**
		 * Tracks the given output file, deleting it if there is an error. 
		 * @param output file to track
		 * @return tracked file
		 */
		public File track(File output) {
			deleteOnError.add(output);
			return output;
		}
	}
}