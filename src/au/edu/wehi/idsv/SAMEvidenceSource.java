package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.EnumSet;
import java.util.Iterator;

import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.metrics.RelevantMetrics;
import au.edu.wehi.idsv.util.AutoClosingIterator;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

/**
 * Structural variation evidence based on read pairs from a single SAM/BAM.  
 * @author cameron.d
 *
 */
public class SAMEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(SAMEvidenceSource.class);
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
	 * @return returns whether any processing was performed
	 */
	public EnumSet<ProcessStep> completeSteps(EnumSet<ProcessStep> steps) {
		ExtractEvidence extract = null;
		try {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)
	    			|| !isComplete(ProcessStep.EXTRACT_SOFT_CLIPS)
	    			|| !isComplete(ProcessStep.EXTRACT_READ_PAIRS)
	    			|| !isComplete(ProcessStep.EXTRACT_READ_MATES)
	    			|| !isComplete(ProcessStep.SORT_READ_MATES)) {
				extract = new ExtractEvidence(processContext, this);
				log.info("START extract evidence for ", input);
				extract.process();
				log.info("SUCCESS extract evidence for ", input);
			}
		}
		catch (Exception e) {
			log.info("FAILURE extract evidence for ", input);
			log.error(e);
		} finally {
			if (extract != null) extract.close();
		}
		return steps;
	}
	@Override
	protected void process() {
		completeSteps(ProcessStep.ALL_STEPS);
	}
	public boolean isComplete(ProcessStep step) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		boolean done = true;
		switch (step) {
			case CALCULATE_METRICS:
				done &= IntermediateFileUtil.checkIntermediate(fsc.getInsertSizeMetrics(input), input);
				done &= IntermediateFileUtil.checkIntermediate(fsc.getIdsvMetrics(input), input);
			case EXTRACT_SOFT_CLIPS:
				if (processContext.shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
						done &= IntermediateFileUtil.checkIntermediate(fsc.getSoftClipBamForChr(input, seq.getSequenceName()), input);
					}
				} else {
					done &= IntermediateFileUtil.checkIntermediate(fsc.getSoftClipBam(input), input);
				}
				break;
			case EXTRACT_READ_PAIRS:
				if (processContext.shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
						done &= IntermediateFileUtil.checkIntermediate(fsc.getReadPairBamForChr(input, seq.getSequenceName()), input);
					}
				} else {
					done &= IntermediateFileUtil.checkIntermediate(fsc.getReadPairBam(input), input);
				}
				break;
			case EXTRACT_READ_MATES:
				if (processContext.shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
						done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBamUnsortedForChr(input, seq.getSequenceName()), input);
					}
				} else {
					done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBamUnsorted(input), input);
				}
				// Unsorted mate files are deleted after SORT_READ_MATES is complete
				done |= isComplete(ProcessStep.SORT_READ_MATES);
				break;
			case SORT_READ_MATES:
				if (processContext.shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
						done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBamForChr(input, seq.getSequenceName()), input);
					}
				} else {
					done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBam(input), input);
				}
				break;
			case REALIGN_SOFT_CLIPS:
				done = isRealignmentComplete();
				break;
			case SORT_REALIGNED_SOFT_CLIPS:
				// FIXME: sort soft clip realignment
				done = false;
			default:
				done = false;
		}
		return done;
	}
	public RelevantMetrics getMetrics() {
		if (metrics == null) {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)) {
				process();
			}
			metrics = new IdsvSamFileMetrics(processContext.getFileSystemContext().getInsertSizeMetrics(input), processContext.getFileSystemContext().getIdsvMetrics(input));
		}
		return metrics;
	}
	public File getSourceFile() { return input; }
	public boolean isTumour() { return isTumour; }
	@Override
	protected Iterator<DirectedEvidence> perChrIterator(String chr) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return iterator(
				fsc.getReadPairBamForChr(input, chr),
				fsc.getMateBamForChr(input, chr),
				fsc.getSoftClipBamForChr(input, chr),
				fsc.getRealignmentBamForChr(input, chr));
	}
	@Override
	protected Iterator<DirectedEvidence> singleFileIterator() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return iterator(
				fsc.getReadPairBam(input),
				fsc.getMateBam(input),
				fsc.getSoftClipBam(input),
				fsc.getRealignmentBam(input));
	}
	private Iterator<DirectedEvidence> iterator(File readPair, File pairMate, File softClip, File realigned) {
		if (!isComplete(ProcessStep.EXTRACT_SOFT_CLIPS) ||
			!isComplete(ProcessStep.EXTRACT_READ_PAIRS) ||
			!isComplete(ProcessStep.EXTRACT_READ_MATES) ||
			!isComplete(ProcessStep.SORT_READ_MATES)) {
			throw new IllegalStateException("Cannot traverse evidence before evidence extraction");
		}
		Iterator<NonReferenceReadPair> rpIt = new ReadPairEvidenceIterator(this,
				processContext.getSamReaderIterator(processContext.getSamReader(readPair)),
				processContext.getSamReaderIterator(processContext.getSamReader(pairMate)));
		rpIt = new AutoClosingIterator<NonReferenceReadPair>(rpIt);
		// sort as input is ordered by SAMRecord alignment start position (SAM/BAM genomic coordinate sorted)
		rpIt = new DirectEvidenceWindowedSortingIterator<NonReferenceReadPair>(processContext, getMetrics().getMaxFragmentSize(), rpIt);
		Iterator<SoftClipEvidence> scIt = new SoftClipEvidenceIterator(processContext, this,
				processContext.getSamReaderIterator(processContext.getSamReader(softClip)));
		scIt = new AutoClosingIterator<SoftClipEvidence>(scIt);
		if (isRealignmentComplete()) {
			scIt = new RealignedSoftClipEvidenceIterator(scIt,
				processContext.getSamReaderIterator(processContext.getSamReader(realigned)));
			scIt = new AutoClosingIterator<SoftClipEvidence>(scIt);
		} else {
			log.info(String.format("Soft clip realignment for %s not completed", softClip));
		}
		// sort into evidence order
		scIt = new DirectEvidenceWindowedSortingIterator<SoftClipEvidence>(processContext, getMetrics().getMaxReadLength(), scIt);
		return Iterators.mergeSorted(ImmutableList.of(rpIt, scIt), DirectedEvidenceOrder.ByNatural);
	}
}