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
import java.util.EnumSet;
import java.util.List;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.FastqBreakpointWriter;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.IdsvMetricsCollector;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.metrics.InsertSizeDistribution;
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
	 * @param steps 
	 */
	private void deleteOutput(EnumSet<ProcessStep> steps) {
		close(); // close any file handles that are still around
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (steps.contains(ProcessStep.CALCULATE_METRICS)) tryDelete(fsc.getInsertSizeMetrics(source.getSourceFile()));
		if (steps.contains(ProcessStep.CALCULATE_METRICS)) tryDelete(fsc.getIdsvMetrics(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getReadPairBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBamUnsorted(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) tryDelete(fsc.getSoftClipBam(source.getSourceFile()));
		if (steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) tryDelete(fsc.getRealignmentFastq(source.getSourceFile()));
		for (SAMSequenceRecord seqr : processContext.getReference().getSequenceDictionary().getSequences()) {
			String seq = seqr.getSequenceName();
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getReadPairBamForChr(source.getSourceFile(), seq));
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBamUnsortedForChr(source.getSourceFile(), seq));
			if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) tryDelete(fsc.getMateBamForChr(source.getSourceFile(), seq));
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
	public void process(final EnumSet<ProcessStep> steps) {
		deleteOutput(steps);
		EnumSet<ProcessStep> firstPassSteps = steps.clone();
		EnumSet<ProcessStep> secondPassSteps = EnumSet.noneOf(ProcessStep.class);
		if (!source.isComplete(ProcessStep.CALCULATE_METRICS)) {
			if (firstPassSteps.remove(ProcessStep.EXTRACT_READ_PAIRS)) {
				secondPassSteps.add(ProcessStep.EXTRACT_READ_PAIRS);
			}
		}
		doProcess(firstPassSteps);
		if (processContext.getReadPairParameters().useProperPairFlag &&
				source.getMetrics().getMaxFragmentSize() > 2 * source.getMetrics().getMedianFragmentSize() + 10 * source.getMetrics().getFragmentSizeStdDev()) {
			log.error(String.format("Proper pair flag indicates fragment size of %d is expected!"
					+ " Use READ_PAIR_CONCORDANT_PERCENT to override discordant pair calculation"
					+ " or realign with an aligner that consider fragment size when setting the proper pair flag.",
					source.getMetrics().getMaxFragmentSize()));
		}
		if (!secondPassSteps.isEmpty()) {
			doProcess(secondPassSteps);
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
	    	
	    	final InsertSizeMetricsCollector ismc = IdsvSamFileMetrics.createInsertSizeMetricsCollector(header);;
	    	final IdsvMetricsCollector imc = IdsvSamFileMetrics.createIdsvMetricsCollector();
	    	
	    	createOutputWriters(steps, header, dictionary);
	    	
	    	processInputRecords(steps, dictionary, ismc, imc);
	    	
	    	flush();
	    	if (steps.contains(ProcessStep.CALCULATE_METRICS)) {
	    		writeMetrics(ismc, imc);
	    	}
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
	private void processInputRecords(final EnumSet<ProcessStep> steps, final SAMSequenceDictionary dictionary, final InsertSizeMetricsCollector ismc, final IdsvMetricsCollector imc) {
		assert((scwriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) || scwriters.size() == dictionary.getSequences().size() || scwriters.size() == 1);
		assert((realignmentWriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_SOFT_CLIPS)) || realignmentWriters.size() == dictionary.getSequences().size() || realignmentWriters.size() == 1);
		assert((rpwriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) || rpwriters.size() == dictionary.getSequences().size() + 1 || rpwriters.size() == 1); // +1 for unmapped 
		assert((matewriters.size() == 0 && !steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) || matewriters.size() == dictionary.getSequences().size() || matewriters.size() == 1);
		
		// Traverse the source.getSourceFile() file
		final ProgressLogger progress = new ProgressLogger(log);
		CloseableIterator<SAMRecord> iter = null;
		try {
			iter = processContext.getSamReaderIterator(reader); 
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
						boolean isConcordant = SAMRecordUtil.couldBeProperPair(record, source.getMetrics().getPairOrientation());
						if (processContext.getReadPairParameters().useProperPairFlag) {
							isConcordant &= record.getProperPairFlag();
						}
						double concordantPercentage = processContext.getReadPairParameters().concordantPercent; // 0.95
						if (concordantPercentage > 0 && isConcordant) {
							InsertSizeDistribution dist = source.getMetrics().getInsertSizeDistribution();
							int insertSize = Math.abs(record.getInferredInsertSize());
							// restrict to reads that are outside the expected distribution
							isConcordant &= dist.cumulativeProbability(insertSize) <= processContext.getReadPairParameters().getCordantPercentageUpperBound()
									&& dist.descendingCumulativeProbability(insertSize) <= processContext.getReadPairParameters().getCordantPercentageUpperBound();
						}
						if (!isConcordant) {
							if (processContext.getReadPairParameters().meetsEvidenceCritera(record)) {
								// don't write this record if we're going to filter
								// still need to write mate as we don't know if it will pass the filter (eg mapq fail here, pass for mate)
								rpwriters.get(offset % rpwriters.size()).addAlignment(record);
							}
							if (!record.getMateUnmappedFlag()) {
								// would be nice to know if we actually need to write this but we don't know enough about the mate yet
								if (steps.contains(ProcessStep.EXTRACT_READ_PAIRS)) matewriters.get(record.getMateReferenceIndex() % matewriters.size()).addAlignment(record);
							}
						}
					}
				}
				if (steps.contains(ProcessStep.CALCULATE_METRICS)) {
					ismc.acceptRecord(record, null);
					imc.acceptRecord(record, null);
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
