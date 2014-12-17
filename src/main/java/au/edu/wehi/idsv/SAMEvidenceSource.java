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
import java.util.Collections;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.pipeline.ExtractEvidence;
import au.edu.wehi.idsv.pipeline.SortRealignedSoftClips;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;

import com.google.common.base.Function;
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
				done &= IntermediateFileUtil.checkIntermediate(fsc.getSoftClipMetrics(input), input);
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
			metrics = new IdsvSamFileMetrics(processContext, getSourceFile());
		}
		return metrics;
	}
	public File getSourceFile() { return input; }
	public boolean isTumour() { return isTumour; }
	public CloseableIterator<DirectedEvidence> iterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote) {
		if (processContext.shouldProcessPerChromosome()) {
			// Lazily iterator over each input
			return new PerChromosomeAggregateIterator<DirectedEvidence>(processContext.getReference().getSequenceDictionary(), new Function<String, Iterator<DirectedEvidence>>() {
				@Override
				public Iterator<DirectedEvidence> apply(String chr) {
					return perChrIterator(includeReadPair, includeSoftClip, includeSoftClipRemote, chr);
				}
			});
		} else {
			return singleFileIterator(includeReadPair, includeSoftClip, includeSoftClipRemote);
		}
	}
	public CloseableIterator<DirectedEvidence> iterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote, final String chr) {
		if (processContext.shouldProcessPerChromosome()) {
			return perChrIterator(includeReadPair, includeSoftClip, includeSoftClipRemote, chr);
		} else {
			return new ChromosomeFilteringIterator<DirectedEvidence>(singleFileIterator(includeReadPair, includeSoftClip, includeSoftClipRemote), processContext.getDictionary().getSequence(chr).getSequenceIndex(), true);
		}
	}
	protected CloseableIterator<DirectedEvidence> perChrIterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote, final String chr) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return iterator(includeReadPair, includeSoftClip, includeSoftClipRemote,
				fsc.getReadPairBamForChr(input, chr),
				fsc.getMateBamForChr(input, chr),
				fsc.getSoftClipBamForChr(input, chr),
				fsc.getRealignmentBamForChr(input, chr),
				fsc.getSoftClipRemoteBamForChr(input, chr),
				fsc.getRealignmentRemoteBamForChr(input, chr),
				chr);
	}
	protected CloseableIterator<DirectedEvidence> singleFileIterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return iterator(includeReadPair, includeSoftClip, includeSoftClipRemote,
				fsc.getReadPairBam(input),
				fsc.getMateBam(input),
				fsc.getSoftClipBam(input),
				fsc.getRealignmentBam(input),
				fsc.getSoftClipRemoteBam(input),
				fsc.getRealignmentRemoteBam(input),
				"");
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
	@SuppressWarnings("unchecked")
	protected CloseableIterator<DirectedEvidence> iterator(
			final boolean includeReadPair,
			final boolean includeSoftClip,
			final boolean includeSoftClipRemote,
			final File readPair,
			final File pairMate,
			final File softClip,
			final File realigned,
			final File remoteSoftClip,
			final File remoteRealigned,
			final String chr) {
		List<CloseableIterator<DirectedEvidence>> itList = Lists.newArrayList();
		if (!isComplete(ProcessStep.EXTRACT_SOFT_CLIPS) ||
			!isComplete(ProcessStep.EXTRACT_READ_PAIRS)) {
			throw new IllegalStateException("Cannot traverse evidence before evidence extraction");
		}
		if (includeReadPair) {
			final CloseableIterator<SAMRecord> rawPairIt = processContext.getSamReaderIterator(readPair);
			final CloseableIterator<SAMRecord> rawMateIt = processContext.getSamReaderIterator(pairMate);
			final CloseableIterator<NonReferenceReadPair> rawRpIt = new ReadPairEvidenceIterator(
					processContext,
					this,
					rawPairIt,
					rawMateIt);
			final Iterator<NonReferenceReadPair> sortedRpIt = new DirectEvidenceWindowedSortingIterator<NonReferenceReadPair>(processContext, getReadPairSortWindowSize(), rawRpIt);
			final CloseableIterator<NonReferenceReadPair> finalrpIt = new AutoClosingIterator<NonReferenceReadPair>(sortedRpIt,
					Lists.<Closeable>newArrayList(rawRpIt, rawPairIt, rawMateIt));
			itList.add((CloseableIterator<DirectedEvidence>)(Object)finalrpIt);
		}
		if (includeSoftClip) {
			final List<Closeable> scToClose = Lists.newArrayList();
			final CloseableIterator<SAMRecord> rawScIt = processContext.getSamReaderIterator(softClip);
			scToClose.add(rawScIt);
			CloseableIterator<SoftClipEvidence> scIt = new SoftClipEvidenceIterator(processContext, this,
					rawScIt);
			scToClose.add(scIt);
			if (isRealignmentComplete()) {
				//log.debug("Realignment is complete for ", input);
				CloseableIterator<SAMRecord> rawRealignIt = processContext.getSamReaderIterator(realigned);
				scToClose.add(rawRealignIt);
				scIt = new RealignedSoftClipEvidenceIterator(scIt, rawRealignIt);
				scToClose.add(scIt);
			} else {
				//log.info("Realigned soft clip evidence not present due to missing realignment bam ", realigned);
			}
			// sort into evidence order
			scIt =  new AutoClosingIterator<SoftClipEvidence>(
						new DirectEvidenceWindowedSortingIterator<SoftClipEvidence>(processContext, getSoftClipSortWindowSize(), scIt),
						scToClose
					);
			itList.add((CloseableIterator<DirectedEvidence>)(Object)scIt);
		}
		if (includeSoftClipRemote) {
			if (!isRealignmentComplete()) {
				throw new IllegalArgumentException(String.format("Realignment resorting not complete for ", input));
			}
			CloseableIterator<RealignedRemoteSoftClipEvidence> remoteScIt = new AutoClosingIterator<RealignedRemoteSoftClipEvidence>(Collections.<RealignedRemoteSoftClipEvidence>emptyIterator()); 
			CloseableIterator<SAMRecord> rsrRawIt = processContext.getSamReaderIterator(remoteRealigned);
			CloseableIterator<SAMRecord> sssRawItf = processContext.getSamReaderIterator(remoteSoftClip);
			CloseableIterator<SAMRecord> sssRawItb = processContext.getSamReaderIterator(remoteSoftClip);
			remoteScIt = new RealignedRemoteSoftClipEvidenceIterator(processContext, this, rsrRawIt, sssRawItf, sssRawItb);
			remoteScIt = new AutoClosingIterator<RealignedRemoteSoftClipEvidence>(new DirectEvidenceWindowedSortingIterator<RealignedRemoteSoftClipEvidence>(processContext, getSoftClipSortWindowSize(), remoteScIt),
					ImmutableList.<Closeable>of(remoteScIt, rsrRawIt, sssRawItf, sssRawItb));
			itList.add((CloseableIterator<DirectedEvidence>)(Object)remoteScIt);
		}
		final CloseableIterator<DirectedEvidence> mergedIt = new AutoClosingMergedIterator<DirectedEvidence>(itList, DirectedEvidenceOrder.ByNatural);
		return new AsyncBufferedIterator<DirectedEvidence>(mergedIt, input.getName() + "-" + chr);
	}
	private int getSoftClipSortWindowSize() {
		return getMaxReadMappedLength() + 1;
	}
	private int getReadPairSortWindowSize() {
		// worst case scenario is a fully mapped forward read followed by a large soft-clipped backward read at max distance
		// MMMMMMMMMMMMMMMMMM>
		// ^--read length --^
		// ^--------max concordant fragment size-----^
		//                   |-----breakend call-----|
		//                   |----------breakend call--------|
		//                                                   <SSSSSSSSSM
		//                   ^--------max concordant fragment size-----^
		// ^ alignment start                                           ^ alignment start
		return getMaxConcordantFragmentSize() + getMaxReadMappedLength() + 1;
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
		return Math.max(getMaxReadMappedLength(), Math.max(getMaxReadLength(), getReadPairConcordanceCalculator().maxConcordantFragmentSize()));
	}
	public int getMaxReadLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_LENGTH;
	}
	public int getMaxReadMappedLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_MAPPED_LENGTH;
	}
}