package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Iterator;

import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.metrics.RelevantMetrics;
import au.edu.wehi.idsv.pipeline.ExtractEvidence;
import au.edu.wehi.idsv.pipeline.SortRealignedSoftClips;
import au.edu.wehi.idsv.util.AutoClosingIterator;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

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
	private SAMFileHeader header;
	public SAMEvidenceSource(ProcessingContext processContext, File file, boolean isTumour) {
		super(processContext, file);
		this.processContext = processContext;
		this.input = file;
		this.isTumour = isTumour;
	}
	public SAMFileHeader getHeader() {
		if (header == null) {
			SamReader reader = null;
			try {
				reader = processContext.getSamReader(input);
				header = reader.getFileHeader();
			} finally {
				CloserUtil.close(reader);
				reader = null;
			}
		}
		return header;
	}
	/**
	 * Ensures that all structural variation evidence has been extracted from the input file 
	 * @return returns steps that could not be performed
	 */
	public EnumSet<ProcessStep> completeSteps(EnumSet<ProcessStep> steps) {
		ExtractEvidence extract = null;
		try {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)
	    			|| !isComplete(ProcessStep.EXTRACT_SOFT_CLIPS)
	    			|| !isComplete(ProcessStep.EXTRACT_READ_PAIRS)) {
				extract = new ExtractEvidence(processContext, this);
				log.info("START extract evidence for ", input);
				extract.process(steps);
				log.info("SUCCESS extract evidence for ", input);
			}
		}
		catch (Exception e) {
			log.info("FAILURE extract evidence for ", input);
			log.error(e);
		} finally {
			if (extract != null) extract.close();
		}
		if (isRealignmentComplete() && steps.contains(ProcessStep.SORT_REALIGNED_SOFT_CLIPS)) {
			SortRealignedSoftClips srsc = new SortRealignedSoftClips(processContext, this);
			if (!srsc.isComplete()) {
				srsc.process(steps);
			}
			srsc.close();
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
				if (!done) return false;
				done &= IntermediateFileUtil.checkIntermediate(fsc.getIdsvMetrics(input), input);
				if (!done) return false;
				break;
			case EXTRACT_SOFT_CLIPS:
				if (processContext.shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
						done &= IntermediateFileUtil.checkIntermediate(fsc.getSoftClipBamForChr(input, seq.getSequenceName()), input);
						if (!done) return false;
					}
				} else {
					done &= IntermediateFileUtil.checkIntermediate(fsc.getSoftClipBam(input), input);
					if (!done) return false;
				}
				break;
			case EXTRACT_READ_PAIRS:
				if (processContext.shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
						done &= IntermediateFileUtil.checkIntermediate(fsc.getReadPairBamForChr(input, seq.getSequenceName()), input);
						if (!done) return false;
						//done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBamUnsortedForChr(input, seq.getSequenceName()), input);
						done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBamForChr(input, seq.getSequenceName()), input);
						if (!done) return false;
					}
				} else {
					done &= IntermediateFileUtil.checkIntermediate(fsc.getReadPairBam(input), input);
					if (!done) return false;
					//done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBamUnsorted(input), input);
					done &= IntermediateFileUtil.checkIntermediate(fsc.getMateBam(input), input);
					if (!done) return false;
				}
				break;
			case REALIGN_SOFT_CLIPS:
				done = isRealignmentComplete();
				break;
			case SORT_REALIGNED_SOFT_CLIPS:
				done = new SortRealignedSoftClips(processContext, this).isComplete();
				break;
			default:
				done = false;
				break;
		}
		return done;
	}
	public RelevantMetrics getMetrics() {
		if (metrics == null) {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)) {
				process();
			}
			metrics = new IdsvSamFileMetrics(processContext, processContext.getFileSystemContext().getInsertSizeMetrics(input), processContext.getFileSystemContext().getIdsvMetrics(input));
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
				fsc.getRealignmentBamForChr(input, chr),
				fsc.getSoftClipRemoteBamForChr(input, chr),
				fsc.getRealignmentRemoteBamForChr(input, chr));
	}
	@Override
	protected Iterator<DirectedEvidence> singleFileIterator() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return iterator(
				fsc.getReadPairBam(input),
				fsc.getMateBam(input),
				fsc.getSoftClipBam(input),
				fsc.getRealignmentBam(input),
				fsc.getSoftClipRemoteBam(input),
				fsc.getRealignmentRemoteBam(input));
	}
	private Iterator<DirectedEvidence> iterator(File readPair, File pairMate, File softClip, File realigned, File remoteSoftClip, File remoteRealigned) {
		if (!isComplete(ProcessStep.EXTRACT_SOFT_CLIPS) ||
			!isComplete(ProcessStep.EXTRACT_READ_PAIRS)) {
			throw new IllegalStateException("Cannot traverse evidence before evidence extraction");
		}
		SamReader rpReader = processContext.getSamReader(readPair);
		SamReader mateReader = processContext.getSamReader(pairMate);
		Iterator<NonReferenceReadPair> rpIt = new ReadPairEvidenceIterator(
				processContext,
				this,
				processContext.getSamReaderIterator(rpReader),
				processContext.getSamReaderIterator(mateReader));
		rpIt = new AutoClosingIterator<NonReferenceReadPair>(rpIt, Lists.<Closeable>newArrayList(rpReader, mateReader));
		// sort as input is ordered by SAMRecord alignment start position (SAM/BAM genomic coordinate sorted)
		rpIt = new DirectEvidenceWindowedSortingIterator<NonReferenceReadPair>(processContext, getMetrics().getMaxFragmentSize(), rpIt);
		SamReader scReader = processContext.getSamReader(softClip);
		Iterator<SoftClipEvidence> scIt = new SoftClipEvidenceIterator(processContext, this,
				processContext.getSamReaderIterator(scReader));
		scIt = new AutoClosingIterator<SoftClipEvidence>(scIt, Lists.<Closeable>newArrayList(scReader));
		if (isRealignmentComplete()) {
			log.debug("Found realignment bam for ", input);
			SamReader scRealignReader = processContext.getSamReader(realigned);
			scIt = new RealignedSoftClipEvidenceIterator(scIt,
				processContext.getSamReaderIterator(scRealignReader));
			scIt = new AutoClosingIterator<SoftClipEvidence>(scIt, Lists.<Closeable>newArrayList(scRealignReader));
		} else {
			log.info("Realigned soft clip evidence not present due to missing realignment bam ", realigned);
		}
		// sort into evidence order
		scIt = new DirectEvidenceWindowedSortingIterator<SoftClipEvidence>(processContext, getMetrics().getMaxReadLength(), scIt);
		
		Iterator<RealignedRemoteSoftClipEvidence> remoteScIt = Collections.emptyIterator(); 
		if (isRealignmentComplete()) {
			if (!isComplete(ProcessStep.SORT_REALIGNED_SOFT_CLIPS)) {
				log.debug("Realignment resorting not complete for ", input);
			} else {
				SamReader realignSortedReader = processContext.getSamReader(remoteRealigned);
				SamReader scRealignSorted1 = processContext.getSamReader(remoteSoftClip);
				SamReader scRealignSorted2 = processContext.getSamReader(remoteSoftClip);
				remoteScIt = new RealignedRemoteSoftClipEvidenceIterator(processContext, this,
						processContext.getSamReaderIterator(realignSortedReader),
						processContext.getSamReaderIterator(scRealignSorted1),
						processContext.getSamReaderIterator(scRealignSorted2));
				remoteScIt = new AutoClosingIterator<RealignedRemoteSoftClipEvidence>(remoteScIt, Lists.<Closeable>newArrayList(realignSortedReader, scRealignSorted1, scRealignSorted2));
				remoteScIt = new DirectEvidenceWindowedSortingIterator<RealignedRemoteSoftClipEvidence>(processContext, getMetrics().getMaxReadLength(), remoteScIt);
			}
		} else {
			log.info("Realigned remote soft clip evidence not present due to missing realignment bam ", realigned);
		}
		return Iterators.mergeSorted(ImmutableList.of(rpIt, scIt, remoteScIt), DirectedEvidenceOrder.ByNatural);
	}
}