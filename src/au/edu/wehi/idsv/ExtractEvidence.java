package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.List;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;
import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.IdsvMetricsCollector;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordMateCoordinateComparator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.Lists;
/**
 * Extracts reads supporting structural variation into intermediate files.
 * By sorting DP and OEA read pairs by the coordinate of the mapped mate,
 * putative directed breakpoints can be assembled by downstream processes
 * in a single pass of the intermediate files.
 * 
 * @author Daniel Cameron
 *
 */
public class ExtractEvidence implements Closeable {
	private static final Log log = Log.getInstance(ExtractEvidence.class);
	public static final String FLAG_NAME = "extracted";
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	public ExtractEvidence(ProcessingContext processContext, SAMEvidenceSource source) {
		this.processContext = processContext;
		this.source = source;
	}
	private SamReader reader = null;
	private ReferenceSequenceFileWalker referenceWalker = null;
	private final List<SAMFileWriter> scwriters = Lists.newArrayList();
	private final List<SAMFileWriter> rpwriters = Lists.newArrayList();
	private final List<SAMFileWriter> matewriters = Lists.newArrayList();
	private final List<FastqBreakpointWriter> realignmentWriters = Lists.newArrayList();
	public void close() {
		CloserUtil.close(reader);
		reader = null;
		CloserUtil.close(referenceWalker);
		referenceWalker = null;
		for (SAMFileWriter w : scwriters) {
			CloserUtil.close(w);
		}
		scwriters.clear();
		for (SAMFileWriter w : rpwriters) {
			CloserUtil.close(w);
		}
		rpwriters.clear();
    	for (SAMFileWriter w : matewriters) {
			CloserUtil.close(w);
		}
    	matewriters.clear();
    	for (FastqBreakpointWriter w : realignmentWriters) {
			CloserUtil.close(w);
		}
    	realignmentWriters.clear();
	}
	/**
	 * Deletes all output files
	 * Deleting output files if there is an error prevents downstream issues
	 * with partially written intermediate files 
	 */
	private void deleteOutput() {
		close(); // close any file handles that are still around
		FileSystemContext fsc = processContext.getFileSystemContext();
		tryDelete(fsc.getFlagFile(source.getSourceFile(), FLAG_NAME));
		tryDelete(fsc.getInsertSizeMetrics(source.getSourceFile()));
		tryDelete(fsc.getIdsvMetrics(source.getSourceFile()));
		tryDelete(fsc.getReadPairBam(source.getSourceFile()));
		tryDelete(fsc.getSoftClipBam(source.getSourceFile()));
		tryDelete(fsc.getMateBam(source.getSourceFile()));
		tryDelete(fsc.getRealignmentFastq(source.getSourceFile()));
		for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
			String seq = seqr.getSequenceName();
			tryDelete(fsc.getReadPairBamForChr(source.getSourceFile(), seq));
			tryDelete(fsc.getSoftClipBamForChr(source.getSourceFile(), seq));
			tryDelete(fsc.getMateBamForChr(source.getSourceFile(), seq));
			tryDelete(fsc.getRealignmentFastqForChr(source.getSourceFile(), seq));
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
	public void process() {
		boolean shouldDelete = true;
    	try {
    		deleteOutput(); // remove any left-over files
    		
	    	reader = processContext.getSamReader(source.getSourceFile());
	    	final SAMFileHeader header = reader.getFileHeader();
	    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
	    	referenceWalker = new ReferenceSequenceFileWalker(processContext.getReferenceFile());
	    	final InsertSizeMetricsCollector ismc = IdsvSamFileMetrics.createInsertSizeMetricsCollector(header);
	    	final IdsvMetricsCollector imc = IdsvSamFileMetrics.createIdsvMetricsCollector();
	    	
	    	dictionary.assertSameDictionary(processContext.getReference().getSequenceDictionary());
	    	
	    	createOutputWriters(header, dictionary);
	    	
	    	processInputRecords(dictionary, ismc, imc);
	    	
	    	flush();
	    	
	    	writeMetrics(ismc, imc);
	    	
	    	sortMates();
	    	
	    	if (!processContext.getFileSystemContext().getFlagFile(source.getSourceFile(), FLAG_NAME).createNewFile()) {
	    		log.error("Failed to create flag file ", processContext.getFileSystemContext().getFlagFile(source.getSourceFile(), FLAG_NAME));
	    	} else {
	    		shouldDelete = false;
	    	}
    	} catch (IOException e) {
    		log.error(e);
			e.printStackTrace();
			shouldDelete = true;
			throw new RuntimeException(e);
		} finally {
    		close();
    		if (shouldDelete) deleteOutput();
    		// Remove temp files
    		tryDelete(processContext.getFileSystemContext().getMateBamUnsorted(source.getSourceFile()));
			for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
				String seq = seqr.getSequenceName();
				tryDelete(processContext.getFileSystemContext().getMateBamUnsortedForChr(source.getSourceFile(), seq));
			}
    	}
    }
	private void sortMates() {
		log.info("Sorting sv mates");
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
				String seq = seqr.getSequenceName();
				SAMFileUtil.sort(processContext,
						processContext.getFileSystemContext().getMateBamUnsortedForChr(source.getSourceFile(), seq),
						processContext.getFileSystemContext().getMateBamForChr(source.getSourceFile(), seq),
						new SAMRecordMateCoordinateComparator());
			}
		} else {
			SAMFileUtil.sort(processContext,
					processContext.getFileSystemContext().getMateBamUnsorted(source.getSourceFile()),
					processContext.getFileSystemContext().getMateBam(source.getSourceFile()),
					new SAMRecordMateCoordinateComparator());
		}
	}
	private void writeMetrics(InsertSizeMetricsCollector ismc, IdsvMetricsCollector imc) {
		log.info("Writing metrics");
		
		ismc.finish();
		MetricsFile<InsertSizeMetrics, Integer> ismmf = processContext.<InsertSizeMetrics, Integer>createMetricsFile();
		ismc.addAllLevelsToFile(ismmf);
		ismmf.write(processContext.getFileSystemContext().getInsertSizeMetrics(source.getSourceFile()));
		
		imc.finish();
		MetricsFile<IdsvMetrics, Integer> imcf = processContext.<IdsvMetrics, Integer>createMetricsFile();
		imc.addAllLevelsToFile(imcf);
		imcf.write(processContext.getFileSystemContext().getIdsvMetrics(source.getSourceFile()));
	}
	private void flush() throws IOException {
		reader.close();
		reader = null;
		for (SAMFileWriter w : scwriters) {
			w.close();
		}
		scwriters.clear();
		for (SAMFileWriter w : rpwriters) {
			w.close();
		}
		rpwriters.clear();
		for (SAMFileWriter w : matewriters) {
			w.close();
		}
		matewriters.clear();
		for (FastqBreakpointWriter w : realignmentWriters) {
			w.close();
		}
		realignmentWriters.clear();
	}
	private void processInputRecords(final SAMSequenceDictionary dictionary, final InsertSizeMetricsCollector ismc, final IdsvMetricsCollector imc) {
		assert(scwriters.size() == dictionary.getSequences().size() || scwriters.size() == 1);
		assert(rpwriters.size() == dictionary.getSequences().size() + 1 || rpwriters.size() == 1); // +1 for unmapped 
		assert(matewriters.size() == dictionary.getSequences().size() || matewriters.size() == 1);
		assert(realignmentWriters.size() == dictionary.getSequences().size() || realignmentWriters.size() == 1);
		
		// Traverse the source.getSourceFile() file
		final ProgressLogger progress = new ProgressLogger(log);
		CloseableIterator<SAMRecord> iter = null;
		try {
			iter = processContext.getSamReaderIterator(reader); 
			while (iter.hasNext()) {
				SAMRecord record = iter.next();
				SAMRecordUtil.ensureNmTag(referenceWalker, record);
				int offset = record.getReadUnmappedFlag() ? dictionary.size() : record.getReferenceIndex();
				SoftClipEvidence startEvidence = null;
				SoftClipEvidence endEvidence = null;
				if (SAMRecordUtil.getStartSoftClipLength(record) > 0) {
					startEvidence = SoftClipEvidence.create(processContext, source, BreakendDirection.Backward, record);
					if (processContext.getSoftClipParameters().meetsEvidenceCritera(startEvidence)) {
						if (processContext.getRealignmentParameters().shouldRealignBreakend(startEvidence)) {
							realignmentWriters.get(offset % realignmentWriters.size()).write(startEvidence);
						}
					} else {
						startEvidence = null;
					}
				}
				if (SAMRecordUtil.getEndSoftClipLength(record) > 0) {
					endEvidence = SoftClipEvidence.create(processContext, source, BreakendDirection.Forward, record);
					if (processContext.getSoftClipParameters().meetsEvidenceCritera(endEvidence)) {
						if (processContext.getRealignmentParameters().shouldRealignBreakend(endEvidence)) {
							realignmentWriters.get(offset % realignmentWriters.size()).write(endEvidence);
						}
					} else {
						endEvidence = null;
					}
				}
				if (startEvidence != null || endEvidence != null) {
					scwriters.get(offset % scwriters.size()).addAlignment(record);
				}
				if (SAMRecordUtil.isPartOfNonReferenceReadPair(record)) {
					rpwriters.get(offset % rpwriters.size()).addAlignment(record);
					if (!record.getMateUnmappedFlag()) {
						matewriters.get(record.getMateReferenceIndex() % matewriters.size()).addAlignment(record);
					}
				}
				ismc.acceptRecord(record, null);
				imc.acceptRecord(record, null);
				progress.record(record);
			}
		} finally {
			if (iter != null) iter.close();
		}
	}
	private void createOutputWriters(final SAMFileHeader header, final SAMSequenceDictionary dictionary) {
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
	private void createOutputWriterPerGenome(final SAMFileHeader svHeader, final SAMFileHeader mateHeader) {
		// all writers map to the same one
		scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSoftClipBam(source.getSourceFile())));
		rpwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBam(source.getSourceFile())));
		matewriters.add(processContext.getSamReaderWriterFactory().makeBAMWriter(mateHeader, true, processContext.getFileSystemContext().getMateBamUnsorted(source.getSourceFile()), 0));
		realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastq(source.getSourceFile()))));
	}
	private void createOutputWritersPerChromosome(
			final SAMSequenceDictionary dictionary,
			final SAMFileHeader svHeader,
			final SAMFileHeader mateHeader) {
		for (SAMSequenceRecord seq : dictionary.getSequences()) {
			scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSoftClipBamForChr(source.getSourceFile(), seq.getSequenceName())));
			rpwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBamForChr(source.getSourceFile(), seq.getSequenceName())));
			matewriters.add(processContext.getSamReaderWriterFactory().makeBAMWriter(mateHeader, true, processContext.getFileSystemContext().getMateBamUnsortedForChr(source.getSourceFile(), seq.getSequenceName()), 0));
			realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastqForChr(source.getSourceFile(), seq.getSequenceName()))));
		}
		rpwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBamForChr(source.getSourceFile(), "unmapped")));
	}
}
