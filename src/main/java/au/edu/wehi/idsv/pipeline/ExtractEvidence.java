package au.edu.wehi.idsv.pipeline;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
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

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.FastqBreakpointWriter;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SequentialCoverageThreshold;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetricsCollector;
import au.edu.wehi.idsv.model.Models;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordCigarCleaningIterator;
import au.edu.wehi.idsv.sam.SAMRecordMateCoordinateComparator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;

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
	private boolean mapqWarned = false;
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
		List<File> toDelete = new ArrayList<>();
		if (steps.contains(ProcessStep.CALCULATE_METRICS)) {
			toDelete.add(fsc.getInsertSizeMetrics(source.getSourceFile()));
			toDelete.add(fsc.getIdsvMetrics(source.getSourceFile()));
			toDelete.add(fsc.getCigarMetrics(source.getSourceFile()));
			toDelete.add(fsc.getMapqMetrics(source.getSourceFile()));
		}
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) toDelete.add(fsc.getReadPairBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) toDelete.add(fsc.getMateBamUnsorted(source.getSourceFile()));
		//if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) toDelete.add(fsc.getSoftClipBam(source.getSourceFile()));
		for (int i = 0; i < source.getRealignmentIterationCount(); i++) {
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) toDelete.add(fsc.getRealignmentFastq(source.getSourceFile(), i));
		}
		for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
			String seq = seqr.getSequenceName();
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) toDelete.add(fsc.getReadPairBamForChr(source.getSourceFile(), seq));
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) toDelete.add(fsc.getMateBamUnsortedForChr(source.getSourceFile(), seq));
			//if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBamForChr(source.getSourceFile(), seq));
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) toDelete.add(fsc.getSoftClipBamForChr(source.getSourceFile(), seq));
			for (int i = 0; i < source.getRealignmentIterationCount(); i++) {
				if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) toDelete.add(fsc.getRealignmentFastqForChr(source.getSourceFile(), seq, i));
			}
		}
		Exception ex = null;
		for (File f : toDelete) {
			if (f.exists()) {
				try {
					if (!f.delete()) {
						log.error("Unable to delete intermediate file ", f,  " during rollback.");
					}
				} catch (Exception e) {
					log.error(e, "Unable to delete intermediate file ", f,  " during rollback.");
					if (ex != null) ex = e;
				}
			}
		}
		if (ex != null) {
			String msg = "Unable to recover from error when rolling back incomplete operation. "
					+ "Invalid intermediate files may exist. "
					+ "Please delete all files in " + fsc.getIntermediateDirectory(source.getSourceFile()).getAbsolutePath();
			log.error(msg);
			throw new RuntimeException(msg, ex);
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
			final IdsvSamFileMetricsCollector collector = new IdsvSamFileMetricsCollector(header);
			final SequentialCoverageThreshold coverageThresholdBlacklist = new SequentialCoverageThreshold(processContext.getDictionary(), processContext.getLinear(), processContext.getConfig().maxCoverage + 1);
	    	
	    	long recordsProcessed = processMetrics(collector, maxRecords, coverageThresholdBlacklist);
	    	writeMetrics(collector);
	    	
	    	if (recordsProcessed >= maxRecords) {
	    		log.info(String.format("Library metrics calculated from first %d records", recordsProcessed));
	    	}
	    	flush();
	    	coverageThresholdBlacklist.finish().write(processContext.getFileSystemContext().getCoverageBlacklistBed(source.getSourceFile()), "coverageBlacklist");
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
	    	if (reader.getFileHeader().getSortOrder() != SortOrder.coordinate) {
	    		String msg = String.format("Input file is not %s coordinate sorted. Please ensure all input files are sorted by read alignment position.", source.getSourceFile());
	    		log.error(msg);
	    		throw new AssertionError(msg);
			}
	    	
	    	// Need to pass the underlying file and reread it since ReferenceSequenceFileWalker forcably
	    	// closes the passed in reference regardless of constructor and we still need it
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
		SAMFileUtil.sort(processContext.getFileSystemContext(), input, output, new SAMRecordMateCoordinateComparator());
		if (deleteInputAfterProcessing) {
			input.delete();
		}
	}
	private void writeMetrics(IdsvSamFileMetricsCollector collector) {
		log.info("Writing metrics");
		collector.finish(processContext, source.getSourceFile());
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
	private long processMetrics(IdsvSamFileMetricsCollector collector, long maxRecords, SequentialCoverageThreshold blacklist) {
		long recordsProcessed = 0;
		final ProgressLogger progress = new ProgressLogger(log, 10000000);
		CloseableIterator<SAMRecord> iter = null;
		try {
			iter = new AsyncBufferedIterator<SAMRecord>(
					new SAMRecordCigarCleaningIterator(processContext.getSamReaderIterator(reader)),
					source.getSourceFile().getName() + "-Metrics",
					processContext.getConfig().async_bufferCount,
					processContext.getConfig().async_bufferSize); 
			while (iter.hasNext() && recordsProcessed++ < maxRecords) {
				SAMRecord record = iter.next();
				collector.acceptRecord(record, null);
				blacklist.acceptRecord(record);
				progress.record(record);
				if (mapqWarned && !record.getReadUnmappedFlag() && record.getMappingQuality() > source.getContext().getConfig().maxMapq) {
					mapqWarned = true;
					log.warn(String.format("Found unexpected mapping quality of %d in %s. Treating all mapping qualities above %d as zero", record.getMappingQuality(), source.getSourceFile(), source.getContext().getConfig().maxMapq));
				}
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
		final ProgressLogger progress = new ProgressLogger(log, 10000000);
		CloseableIterator<SAMRecord> iter = null;
		try {
			iter = new AsyncBufferedIterator<SAMRecord>(
					new SAMRecordCigarCleaningIterator(processContext.getSamReaderIterator(reader)),
					source.getSourceFile().getName() + "-Extract",
					processContext.getConfig().async_bufferCount,
					processContext.getConfig().async_bufferSize);
			while (iter.hasNext()) {
				SAMRecord record = iter.next();
				SAMRecordUtil.ensureNmTag(referenceWalker, record);
				int offset = record.getReadUnmappedFlag() ? dictionary.size() : record.getReferenceIndex();
				if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) {
					if (Models.meetsReadAlignmentCriteria(processContext.getConfig(), record)) {
						SoftClipEvidence startEvidence = null;
						SoftClipEvidence endEvidence = null;
						if (SAMRecordUtil.getStartSoftClipLength(record) > 0) {
							startEvidence = SoftClipEvidence.create(source, BreakendDirection.Backward, record);
							if (startEvidence.meetsEvidenceCritera()) {
								if (processContext.getRealignmentParameters().shouldRealignBreakend(startEvidence)) {
									realignmentWriters.get(offset % realignmentWriters.size()).write(startEvidence);
								}
							} else {
								startEvidence = null;
							}
						}
						if (SAMRecordUtil.getEndSoftClipLength(record) > 0) {
							endEvidence = SoftClipEvidence.create(source, BreakendDirection.Forward, record);
							if (endEvidence.meetsEvidenceCritera()) {
								if (processContext.getRealignmentParameters().shouldRealignBreakend(endEvidence)) {
									realignmentWriters.get(offset % realignmentWriters.size()).write(endEvidence);
								}
							} else {
								endEvidence = null;
							}
						}
						if (startEvidence != null || endEvidence != null || CigarUtil.countIndels(record.getCigar()) > 0) {
							scwriters.get(offset % scwriters.size()).addAlignment(record);
						}
					}
				}
				if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) {
					if (record.getReadPairedFlag()) {
						if (NonReferenceReadPair.meetsAnchorCriteria(source, record)) {
							rpwriters.get(offset % rpwriters.size()).addAlignment(record);
						}
						if (NonReferenceReadPair.meetsRemoteCriteria(source, record)) {
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
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) scwriters.add(processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSoftClipBam(source.getSourceFile())));
		for (int i = 0; i < source.getRealignmentIterationCount(); i++) {
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastq(source.getSourceFile(), i))));
		}
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) rpwriters.add(processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBam(source.getSourceFile())));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) matewriters.add(processContext.getSamFileWriterFactory(false).makeBAMWriter(mateHeader, true, processContext.getFileSystemContext().getMateBamUnsorted(source.getSourceFile()), 0));
		
	}
	private void createOutputWritersPerChromosome(
			final EnumSet<ProcessStep> steps,
			final SAMSequenceDictionary dictionary,
			final SAMFileHeader svHeader,
			final SAMFileHeader mateHeader) {
		for (SAMSequenceRecord seq : dictionary.getSequences()) {
			if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) scwriters.add(processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getSoftClipBamForChr(source.getSourceFile(), seq.getSequenceName())));
			for (int i = 0; i < source.getRealignmentIterationCount(); i++) {
				if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) realignmentWriters.add(new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(processContext.getFileSystemContext().getRealignmentFastqForChr(source.getSourceFile(), seq.getSequenceName(), i))));
			}
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) rpwriters.add(processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBamForChr(source.getSourceFile(), seq.getSequenceName())));
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) matewriters.add(processContext.getSamFileWriterFactory(false).makeBAMWriter(mateHeader, true, processContext.getFileSystemContext().getMateBamUnsortedForChr(source.getSourceFile(), seq.getSequenceName()), 0));
			
		}
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) rpwriters.add(processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(svHeader, true, processContext.getFileSystemContext().getReadPairBamForChr(source.getSourceFile(), "unmapped")));
	}
}
