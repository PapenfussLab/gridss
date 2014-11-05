package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.pipeline.ExtractEvidence;
import au.edu.wehi.idsv.pipeline.SortRealignedSoftClips;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;

import com.google.common.collect.ImmutableList;
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
	private final ReadPairConcordanceMethod rpcMethod;
	private final Object[] rpcParams;
	private IdsvSamFileMetrics metrics;
	private SAMFileHeader header;
	private ReadPairConcordanceCalculator rpcc;
	public SAMEvidenceSource(ProcessingContext processContext, File file, boolean isTumour) {
		this(processContext, file, isTumour, ReadPairConcordanceMethod.SAM_FLAG, null);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, boolean isTumour, int minFragmentSize, int maxFragmentSize) {
		this(processContext, file, isTumour, ReadPairConcordanceMethod.FIXED, new Object[] { minFragmentSize, maxFragmentSize });
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, boolean isTumour, double concordantPercentage) {
		this(processContext, file, isTumour, ReadPairConcordanceMethod.PERCENTAGE, new Object[] { (double)concordantPercentage });
	}
	protected SAMEvidenceSource(ProcessingContext processContext, File file, boolean isTumour, ReadPairConcordanceMethod rpcMethod, Object[] rpcParams) {
		super(processContext, file);
		this.processContext = processContext;
		this.input = file;
		this.isTumour = isTumour;
		this.rpcMethod = rpcMethod;
		this.rpcParams = rpcParams;
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
	public void completeSteps(EnumSet<ProcessStep> steps) {
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
	protected IdsvSamFileMetrics getMetrics() {
		if (metrics == null) {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)) {
				completeSteps(EnumSet.of(ProcessStep.CALCULATE_METRICS));
			}
			metrics = new IdsvSamFileMetrics(processContext.getFileSystemContext().getInsertSizeMetrics(input), processContext.getFileSystemContext().getIdsvMetrics(input));
		}
		return metrics;
	}
	public File getSourceFile() { return input; }
	public boolean isTumour() { return isTumour; }
	@Override
	protected CloseableIterator<DirectedEvidence> perChrIterator(String chr) {
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
	protected CloseableIterator<DirectedEvidence> singleFileIterator() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return iterator(
				fsc.getReadPairBam(input),
				fsc.getMateBam(input),
				fsc.getSoftClipBam(input),
				fsc.getRealignmentBam(input),
				fsc.getSoftClipRemoteBam(input),
				fsc.getRealignmentRemoteBam(input));
	}
	/**
	 * Iterates over the given evidence in evidence breakend position order 
	 * @param readPair
	 * @param pairMate
	 * @param softClip
	 * @param realigned
	 * @param remoteSoftClip
	 * @param remoteRealigned
	 * @return
	 */
	private CloseableIterator<DirectedEvidence> iterator(File readPair, File pairMate, File softClip, File realigned, File remoteSoftClip, File remoteRealigned) {
		if (!isComplete(ProcessStep.EXTRACT_SOFT_CLIPS) ||
			!isComplete(ProcessStep.EXTRACT_READ_PAIRS)) {
			throw new IllegalStateException("Cannot traverse evidence before evidence extraction");
		}
		SamReader rpReader = processContext.getSamReader(readPair);
		SamReader mateReader = processContext.getSamReader(pairMate);
		CloseableIterator<SAMRecord> rawPairIt = processContext.getSamReaderIterator(rpReader);
		CloseableIterator<SAMRecord> rawMateIt = processContext.getSamReaderIterator(mateReader);
		final CloseableIterator<NonReferenceReadPair> rawRpIt = new ReadPairEvidenceIterator(
				processContext,
				this,
				rawPairIt,
				rawMateIt);
		final Iterator<NonReferenceReadPair> sortedRpIt = new DirectEvidenceWindowedSortingIterator<NonReferenceReadPair>(processContext, getMaxConcordantFragmentSize(), rawRpIt);
		final CloseableIterator<NonReferenceReadPair> finalrpIt = new AutoClosingIterator<NonReferenceReadPair>(sortedRpIt,
				Lists.<Closeable>newArrayList(rawRpIt, rpReader, mateReader, rawPairIt, rawMateIt));
		
		List<Closeable> scToClose = new ArrayList<>();
		final SamReader scReader = processContext.getSamReader(softClip);
		scToClose.add(scReader);
		CloseableIterator<SAMRecord> rawScIt = processContext.getSamReaderIterator(scReader);
		scToClose.add(rawScIt);
		CloseableIterator<SoftClipEvidence> scIt = new SoftClipEvidenceIterator(processContext, this,
				rawScIt);
		scToClose.add(scIt);
		if (isRealignmentComplete()) {
			//log.debug("Realignment is complete for ", input);
			SamReader scRealignReader = processContext.getSamReader(realigned);
			scToClose.add(scRealignReader);
			CloseableIterator<SAMRecord> rawRealignIt = processContext.getSamReaderIterator(scRealignReader);
			scToClose.add(rawRealignIt);
			scIt = new RealignedSoftClipEvidenceIterator(scIt,
					rawRealignIt);
			scToClose.add(scIt);
		} else {
			//log.info("Realigned soft clip evidence not present due to missing realignment bam ", realigned);
		}
		// sort into evidence order
		scIt =  new AutoClosingIterator<SoftClipEvidence>(
					new DirectEvidenceWindowedSortingIterator<SoftClipEvidence>(processContext, getMetrics().getIdsvMetrics().MAX_READ_LENGTH, scIt),
					scToClose
				);
		
		CloseableIterator<RealignedRemoteSoftClipEvidence> remoteScIt = new AutoClosingIterator<RealignedRemoteSoftClipEvidence>(Collections.<RealignedRemoteSoftClipEvidence>emptyIterator()); 
		if (isRealignmentComplete()) {
			if (!isComplete(ProcessStep.SORT_REALIGNED_SOFT_CLIPS)) {
				//log.debug("Realignment resorting not complete for ", input);
			} else {
				SamReader realignSortedReader = processContext.getSamReader(remoteRealigned);
				SamReader scRealignSortedf = processContext.getSamReader(remoteSoftClip);
				SamReader scRealignSortedb = processContext.getSamReader(remoteSoftClip);
				CloseableIterator<SAMRecord> rsrRawIt = processContext.getSamReaderIterator(realignSortedReader);
				CloseableIterator<SAMRecord> sssRawItf = processContext.getSamReaderIterator(scRealignSortedf);
				CloseableIterator<SAMRecord> sssRawItb = processContext.getSamReaderIterator(scRealignSortedb);
				remoteScIt = new RealignedRemoteSoftClipEvidenceIterator(processContext, this,
						rsrRawIt, sssRawItf, sssRawItb);
				remoteScIt = new AutoClosingIterator<RealignedRemoteSoftClipEvidence>(new DirectEvidenceWindowedSortingIterator<RealignedRemoteSoftClipEvidence>(processContext, getMetrics().getIdsvMetrics().MAX_READ_LENGTH, remoteScIt),
						ImmutableList.<Closeable>of(remoteScIt, realignSortedReader, scRealignSortedf, scRealignSortedb, rsrRawIt, sssRawItf, sssRawItb));
			}
		} else {
			//log.info("Realigned remote soft clip evidence not present due to missing realignment bam ", realigned);
		}
		return new AutoClosingMergedIterator<DirectedEvidence>(ImmutableList.of(finalrpIt, scIt, remoteScIt), DirectedEvidenceOrder.ByNatural);
	}
	public ReadPairConcordanceCalculator getReadPairConcordanceCalculator() {
		if (rpcc == null) {
			switch (rpcMethod) {
				case FIXED:
					rpcc = new FixedSizeReadPairConcordanceCalculator((int)(Integer)rpcParams[0], (int)(Integer)rpcParams[1]);
					break;
				case PERCENTAGE:
					rpcc = new PercentageReadPairConcordanceCalculator(getMetrics().getInsertSizeDistribution(), (double)(Double)rpcParams[0]);
					break;
				default:
				case SAM_FLAG:
					rpcc = new SAMFlagReadPairConcordanceCalculator(getMetrics().getIdsvMetrics());
					// Safety check for BWA which sets proper pair flag based only correct chromosome and orientation 
					if (getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH >= 100000 &&  
							getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH >= getMetrics().getInsertSizeMetrics().MAX_INSERT_SIZE / 2 &&
							getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH >= 10 * getMaxReadLength()) {
						String msg = String.format("Proper pair flag indicates fragment size of %d is expected!"
								+ " Realign with an aligner that consider fragment size when setting the proper pair flag or "
								+ " specify fixed or percentage bounds for read pair concordance for %s."
								+ " Insert size distribution counts can be found in %s",
								getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH,
								input,
								processContext.getFileSystemContext().getIdsvMetrics(input)); 
						log.error(msg);
						throw new IllegalArgumentException(msg);
					}
					break;
			}
		}
		return rpcc;
	}
	@Override
	public int getMaxConcordantFragmentSize() {
		return Math.max(getMetrics().getIdsvMetrics().MAX_READ_LENGTH, getReadPairConcordanceCalculator().maxConcordantFragmentSize());
	}
	public int getMaxReadLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_LENGTH;
	}
}