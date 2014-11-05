package au.edu.wehi.idsv.pipeline;

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
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.FastqBreakpointWriter;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.IdsvMetricsCollector;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordMateCoordinateComparator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;

import com.google.common.io.Files;
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
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	private final FileSystemContext fsc;
	public ExtractEvidence(ProcessingContext processContext, SAMEvidenceSource source) {
		this.processContext = processContext;
		this.source = source;
		this.fsc = processContext.getFileSystemContext();
	}
	private SamReader reader = null;
	private ReferenceSequenceFileWalker referenceWalker = null;
	private final List<SAMFileWriter> scwriters = new ArrayList<>();
	private final List<SAMFileWriter> rpwriters = new ArrayList<>();
	private final List<SAMFileWriter> matewriters = new ArrayList<>();
	private final List<FastqBreakpointWriter> realignmentWriters = new ArrayList<>();
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
	 * @param steps 
	 */
	private void deleteOutput(EnumSet<ProcessStep> steps) {
		close(); // close any file handles that are still around
		FileSystemContext fsc = processContext.getFileSystemContext();
		//if (steps.contains(ProcessStep.CALCULATE_METRICS)) tryDelete(fsc.getInsertSizeMetrics(source.getSourceFile()));
		//if (steps.contains(ProcessStep.CALCULATE_METRICS)) tryDelete(fsc.getIdsvMetrics(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getReadPairBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBamUnsorted(source.getSourceFile()));
		//if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) tryDelete(fsc.getSoftClipBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) tryDelete(fsc.getRealignmentFastq(source.getSourceFile()));
		for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
			String seq = seqr.getSequenceName();
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getReadPairBamForChr(source.getSourceFile(), seq));
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBamUnsortedForChr(source.getSourceFile(), seq));
			//if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBamForChr(source.getSourceFile(), seq));
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) tryDelete(fsc.getSoftClipBamForChr(source.getSourceFile(), seq));
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) tryDelete(fsc.getRealignmentFastqForChr(source.getSourceFile(), seq));
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
	public void process(EnumSet<ProcessStep> steps) {
		deleteOutput(steps);
		EnumSet<ProcessStep> remainingSteps = steps.clone();
		if (!source.isComplete(ProcessStep.CALCULATE_METRICS) && steps.contains(ProcessStep.CALCULATE_METRICS)) {
			gatherMetrics(processContext.getCalculateMetricsRecordCount());
			remainingSteps.remove(ProcessStep.CALCULATE_METRICS);
		}
		if (!remainingSteps.isEmpty()) {
			doProcess(remainingSteps);
		}
	}
	private void gatherMetrics(long maxRecords) {
		boolean shouldDelete = true;
		try {
			reader = processContext.getSamReader(source.getSourceFile());
			final SAMFileHeader header = reader.getFileHeader();
			final InsertSizeMetricsCollector ismc = IdsvSamFileMetrics.createInsertSizeMetricsCollector(header);
	    	final IdsvMetricsCollector imc = IdsvSamFileMetrics.createIdsvMetricsCollector();
	    	
	    	long recordsProcessed = processMetrics(ismc, imc, maxRecords);
	    	writeMetrics(ismc, imc);
	    	
	    	if (recordsProcessed >= maxRecords) {
	    		log.info(String.format("Library metrics calculated from first %d records", recordsProcessed));
	    	}
	    	flush();
	    	shouldDelete = false;
    	} catch (IOException e) {
    		log.error(e);
			shouldDelete = true;
			throw new RuntimeException(e);
		} finally {
    		close();
    		if (shouldDelete) deleteOutput(EnumSet.of(ProcessStep.CALCULATE_METRICS));
    	}
		
	}
	private void doProcess(EnumSet<ProcessStep> steps) {
		boolean shouldDelete = true;
    	try {
	    	reader = processContext.getSamReader(source.getSourceFile());
	    	final SAMFileHeader header = reader.getFileHeader();
	    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
	    	dictionary.assertSameDictionary(processContext.getReference().getSequenceDictionary());
	    	
	    	referenceWalker = new ReferenceSequenceFileWalker(processContext.getReferenceFile());
	    	
	    	createOutputWriters(steps, header, dictionary);
	    	
	    	processInputRecords(steps, dictionary);
	    	
	    	flush();
	    	if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) {
	    		sortMates();
	    	}
	    	shouldDelete = false;
    	} catch (IOException e) {
    		log.error(e);
			e.printStackTrace();
			shouldDelete = true;
			throw new RuntimeException(e);
		} finally {
    		close();
    		if (shouldDelete) deleteOutput(steps);
    	}
	}
	private void sortMates() {
		log.info("Sorting sv mates");
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
				String seq = seqr.getSequenceName();
				// PARALLEL opportunity - not great candidate due memory usage of and intermediate storage
				// required for parallel sorts
				sortMates(
						processContext.getFileSystemContext().getMateBamUnsortedForChr(source.getSourceFile(), seq),
						processContext.getFileSystemContext().getMateBamForChr(source.getSourceFile(), seq),
						true);
			}
		} else {  
			sortMates(
					processContext.getFileSystemContext().getMateBamUnsorted(source.getSourceFile()),
					processContext.getFileSystemContext().getMateBam(source.getSourceFile()),
					true);
		}
	}
	private void sortMates(File input, File output, boolean deleteInputAfterProcessing) {
		SAMFileUtil.sort(processContext, input, output, new SAMRecordMateCoordinateComparator());
		if (deleteInputAfterProcessing) {
			input.delete();
		}
	}
	private void writeMetrics(InsertSizeMetricsCollector ismc, IdsvMetricsCollector imc) throws IOException {
		log.info("Writing metrics");
		
		ismc.finish();
		MetricsFile<InsertSizeMetrics, Integer> ismmf = processContext.<InsertSizeMetrics, Integer>createMetricsFile();
		ismc.addAllLevelsToFile(ismmf);
		File ismOutput = fsc.getInsertSizeMetrics(source.getSourceFile());
		ismOutput.delete();
		ismmf.write(FileSystemContext.getWorkingFileFor(ismOutput));
		Files.move(FileSystemContext.getWorkingFileFor(ismOutput), ismOutput);
		
		imc.finish();
		MetricsFile<IdsvMetrics, Integer> imcf = processContext.<IdsvMetrics, Integer>createMetricsFile();
		imc.addAllLevelsToFile(imcf);
		File imOutput = fsc.getIdsvMetrics(source.getSourceFile());
		imcf.write(FileSystemContext.getWorkingFileFor(imOutput));
		Files.move(FileSystemContext.getWorkingFileFor(imOutput), imOutput);
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
	private long processMetrics(InsertSizeMetricsCollector ismc, IdsvMetricsCollector imc, long maxRecords) {
		long recordsProcessed = 0;
		final ProgressLogger progress = new ProgressLogger(log);
		CloseableIterator<SAMRecord> iter = null;
		try {
			iter = new AsyncBufferedIterator<SAMRecord>(processContext.getSamReaderIterator(reader), source.getSourceFile().getName() + " metrics"); 
			while (iter.hasNext() && recordsProcessed++ < maxRecords) {
				SAMRecord record = iter.next();
				ismc.acceptRecord(record, null);
				imc.acceptRecord(record, null);
				progress.record(record);
			}
		} finally {
			CloserUtil.close(iter);
		}
		return recordsProcessed;
	}
	private void processInputRecords(final EnumSet<ProcessStep> steps, final SAMSequenceDictionary dictionary) {
		assert((scwriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) || scwriters.size() == dictionary.getSequences().size() || scwriters.size() == 1);
		assert((realignmentWriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) || realignmentWriters.size() == dictionary.getSequences().size() || realignmentWriters.size() == 1);
		assert((rpwriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) || rpwriters.size() == dictionary.getSequences().size() + 1 || rpwriters.size() == 1); // +1 for unmapped 
		assert((matewriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) || matewriters.size() == dictionary.getSequences().size() || matewriters.size() == 1);
		
		// Traverse the source.getSourceFile() file
		final ProgressLogger progress = new ProgressLogger(log);
		CloseableIterator<SAMRecord> iter = null;
		try {
			iter = new AsyncBufferedIterator<SAMRecord>(processContext.getSamReaderIterator(reader), source.getSourceFile().getName() + " extract");
			while (iter.hasNext()) {
				SAMRecord record = iter.next();
				SAMRecordUtil.ensureNmTag(referenceWalker, record);
				int offset = record.getReadUnmappedFlag() ? dictionary.size() : record.getReferenceIndex();
				if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) {
					SoftClipEvidence startEvidence = null;
					SoftClipEvidence endEvidence = null;
					if (SAMRecordUtil.getStartSoftClipLength(record) > 0) {
						startEvidence = SoftClipEvidence.create(processContext, source, BreakendDirection.Backward, record);
						if (startEvidence.meetsEvidenceCritera(processContext.getSoftClipParameters())) {
							if (processContext.getRealignmentParameters().shouldRealignBreakend(startEvidence)) {
								realignmentWriters.get(offset % realignmentWriters.size()).write(startEvidence);
							}
						} else {
							startEvidence = null;
						}
					}
					if (SAMRecordUtil.getEndSoftClipLength(record) > 0) {
						endEvidence = SoftClipEvidence.create(processContext, source, BreakendDirection.Forward, record);
						if (endEvidence.meetsEvidenceCritera(processContext.getSoftClipParameters())) {
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
				}
				if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) {
					if (record.getReadPairedFlag()) {
						if (NonReferenceReadPair.meetsLocalEvidenceCritera(processContext.getReadPairParameters(), source, record)) {
							rpwriters.get(offset % rpwriters.size()).addAlignment(record);
						}
						if (NonReferenceReadPair.meetsRemoteEvidenceCritera(processContext.getReadPairParameters(), source, record)) {
							matewriters.get(record.getMateReferenceIndex() % matewriters.size()).addAlignment(record);
						}
					}
				}
				progress.record(record);
			}
		} finally {
			if (iter != null) iter.close();
		}
	}
	private void createOutputWriters(final EnumSet<ProcessStep> steps, final SAMFileHeader header, final SAMSequenceDictionary dictionary) {
		final SAMFileHeader svHeader = header.clone();
		svHeader.setSortOrder(SortOrder.coordinate);
		final SAMFileHeader mateHeader = header.clone();
		mateHeader.setSortOrder(SortOrder.unsorted);
		if (processContext.shouldProcessPerChromosome()) {
			createOutputWritersPerChromosome(steps, dictionary, svHeader, mateHeader);
		} else {
			createOutputWriterPerGenome(steps, svHeader, mateHeader);
		}
	}
	private void createOutputWriterPerGenome(final EnumSet<ProcessStep> steps, final SAMFileHeader svHeader, final SAMFileHeader mateHeader) {
		// all writers map to the same one
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSoftClipBam(source.getSourceFile())));
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastq(source.getSourceFile()))));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) rpwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBam(source.getSourceFile())));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) matewriters.add(processContext.getSamReaderWriterFactory().makeBAMWriter(mateHeader, true, processContext.getFileSystemContext().getMateBamUnsorted(source.getSourceFile()), 0));
		
	}
	private void createOutputWritersPerChromosome(
			final EnumSet<ProcessStep> steps,
			final SAMSequenceDictionary dictionary,
			final SAMFileHeader svHeader,
			final SAMFileHeader mateHeader) {
		for (SAMSequenceRecord seq : dictionary.getSequences()) {
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) scwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSoftClipBamForChr(source.getSourceFile(), seq.getSequenceName())));
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastqForChr(source.getSourceFile(), seq.getSequenceName()))));
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) rpwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBamForChr(source.getSourceFile(), seq.getSequenceName())));
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) matewriters.add(processContext.getSamReaderWriterFactory().makeBAMWriter(mateHeader, true, processContext.getFileSystemContext().getMateBamUnsortedForChr(source.getSourceFile(), seq.getSequenceName()), 0));
			
		}
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) rpwriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBamForChr(source.getSourceFile(), "unmapped")));
	}
}
