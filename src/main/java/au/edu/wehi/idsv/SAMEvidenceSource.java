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
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.pipeline.ExtractEvidence;
import au.edu.wehi.idsv.pipeline.SortRealignedSoftClips;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

/**
 * Structural variation evidence based on read pairs from a single SAM/BAM.  
 * @author Daniel Cameron
 *
 */
public class SAMEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(SAMEvidenceSource.class);
	private final File input;
	private final int sourceCategory;
	private final ReadPairConcordanceMethod rpcMethod;
	private final Object[] rpcParams;
	private IdsvSamFileMetrics metrics;
	private IntervalBed blacklist;
	private SAMFileHeader header;
	private ReadPairConcordanceCalculator rpcc;
	public SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory) {
		this(processContext, file, sourceCategory, ReadPairConcordanceMethod.SAM_FLAG, null);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory, int minFragmentSize, int maxFragmentSize) {
		this(processContext, file, sourceCategory, ReadPairConcordanceMethod.FIXED, new Object[] { minFragmentSize, maxFragmentSize });
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory, double concordantPercentage) {
		this(processContext, file, sourceCategory, ReadPairConcordanceMethod.PERCENTAGE, new Object[] { (double)concordantPercentage });
	}
	protected SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory, ReadPairConcordanceMethod rpcMethod, Object[] rpcParams) {
		super(processContext, file);
		this.input = file;
		this.sourceCategory = sourceCategory;
		this.rpcMethod = rpcMethod;
		this.rpcParams = rpcParams;
	}
	@Override
	public int getRealignmentIterationCount() {
		return 1;
	}
	public SAMFileHeader getHeader() {
		if (header == null) {
			SamReader reader = null;
			try {
				reader = getContext().getSamReader(input);
				header = reader.getFileHeader();
			} finally {
				CloserUtil.close(reader);
			}
		}
		return header;
	}
	/**
	 * Ensures that all structural variation evidence has been extracted from the input file 
	 */
	public void completeSteps(EnumSet<ProcessStep> steps) {
		ExtractEvidence extract = null;
		try {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)
	    			|| !isComplete(ProcessStep.EXTRACT_SOFT_CLIPS)
	    			|| !isComplete(ProcessStep.EXTRACT_READ_PAIRS)) {
				extract = new ExtractEvidence(getContext(), this);
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
		if (steps.contains(ProcessStep.REALIGN_SOFT_CLIPS)) {
			isRealignmentComplete(true);
		}
		if (isRealignmentComplete(false) && steps.contains(ProcessStep.SORT_REALIGNED_SOFT_CLIPS)) {
			SortRealignedSoftClips srsc = new SortRealignedSoftClips(getContext(), this);
			if (!srsc.isComplete()) {
				srsc.process(steps);
			}
			srsc.close();
		}
	}
	public boolean isComplete(ProcessStep step) {
		FileSystemContext fsc = getContext().getFileSystemContext();
		List<File> target = new ArrayList<File>();
		List<File> source = new ArrayList<File>();
		boolean done = true;
		switch (step) {
			case CALCULATE_METRICS:
				source.add(input); target.add(fsc.getInsertSizeMetrics(input));
				source.add(input); target.add(fsc.getIdsvMetrics(input));
				source.add(input); target.add(fsc.getCigarMetrics(input));
				source.add(input); target.add(fsc.getCoverageBlacklistBed(input));
				break;
			case EXTRACT_SOFT_CLIPS:
				if (getContext().shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : getContext().getReference().getSequenceDictionary().getSequences()) {
						source.add(input); target.add(fsc.getSoftClipBamForChr(input, seq.getSequenceName()));
					}
				} else {
					source.add(input); target.add(fsc.getSoftClipBam(input));
				}
				break;
			case EXTRACT_READ_PAIRS:
				if (getContext().shouldProcessPerChromosome()) {
					for (SAMSequenceRecord seq : getContext().getReference().getSequenceDictionary().getSequences()) {
						source.add(input); target.add(fsc.getReadPairBamForChr(input, seq.getSequenceName()));
						source.add(input); target.add(fsc.getMateBamForChr(input, seq.getSequenceName()));
					}
				} else {
					source.add(input); target.add(fsc.getReadPairBam(input));
					source.add(input); target.add(fsc.getMateBam(input));
				}
				break;
			case REALIGN_SOFT_CLIPS:
				done = isRealignmentComplete(false);
				break;
			case SORT_REALIGNED_SOFT_CLIPS:
				done = new SortRealignedSoftClips(getContext(), this).isComplete();
				break;
			default:
				done = false;
				break;
		}
		return done && IntermediateFileUtil.checkIntermediate(target, source, getContext().getConfig().ignoreFileTimestamps);
	}
	public IdsvSamFileMetrics getMetrics() {
		if (metrics == null) {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)) {
				completeSteps(EnumSet.of(ProcessStep.CALCULATE_METRICS));
			}
			metrics = new IdsvSamFileMetrics(getContext(), getSourceFile());
		}
		return metrics;
	}
	public File getSourceFile() { return input; }
	/**
	 * Grouping category of source. For tumour/normal analysis, normal is considered group 0, tumour group 1
	 * @return
	 */
	public int getSourceCategory() {
		return sourceCategory;
	}
	public CloseableIterator<DirectedEvidence> iterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote) {
		if (getContext().shouldProcessPerChromosome()) {
			// Lazily iterator over each input
			return new PerChromosomeAggregateIterator<DirectedEvidence>(getContext().getReference().getSequenceDictionary(), new Function<String, Iterator<DirectedEvidence>>() {
				@Override
				public Iterator<DirectedEvidence> apply(String chr) {
					return perChrIterator(includeReadPair, includeSoftClip, includeSoftClipRemote, chr);
				}
			});
		} else {
			CloseableIterator<DirectedEvidence> it = singleFileIterator(includeReadPair, includeSoftClip, includeSoftClipRemote);
			// can only check pairing if we are iterating over the entire record set
			if (Defaults.SANITY_CHECK_ITERATORS) {
				if (includeSoftClip == includeSoftClipRemote) {
					it = new AutoClosingIterator<DirectedEvidence>(new PairedEvidenceTracker<DirectedEvidence>(it), ImmutableList.<Closeable>of(it));
				}
			}
			return it; 
		}
	}
	public CloseableIterator<DirectedEvidence> iterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote, final String chr) {
		if (getContext().shouldProcessPerChromosome()) {
			return perChrIterator(includeReadPair, includeSoftClip, includeSoftClipRemote, chr);
		} else {
			return new ChromosomeFilteringIterator<DirectedEvidence>(singleFileIterator(includeReadPair, includeSoftClip, includeSoftClipRemote), getContext().getDictionary().getSequence(chr).getSequenceIndex(), true);
		}
	}
	protected CloseableIterator<DirectedEvidence> perChrIterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote, final String chr) {
		FileSystemContext fsc = getContext().getFileSystemContext();
		return iterator(includeReadPair, includeSoftClip, includeSoftClipRemote,
				fsc.getReadPairBamForChr(input, chr),
				fsc.getMateBamForChr(input, chr),
				fsc.getSoftClipBamForChr(input, chr),
				fsc.getRealignmentBamForChr(input, chr, 0),
				fsc.getSoftClipRemoteBamForChr(input, chr),
				fsc.getRealignmentRemoteBamForChr(input, chr),
				chr);
	}
	protected CloseableIterator<DirectedEvidence> singleFileIterator(final boolean includeReadPair, final boolean includeSoftClip, final boolean includeSoftClipRemote) {
		FileSystemContext fsc = getContext().getFileSystemContext();
		return iterator(includeReadPair, includeSoftClip, includeSoftClipRemote,
				fsc.getReadPairBam(input),
				fsc.getMateBam(input),
				fsc.getSoftClipBam(input),
				fsc.getRealignmentBam(input, 0),
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
	@SuppressWarnings({ "unchecked" })
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
		if (!isComplete(ProcessStep.CALCULATE_METRICS) ||
			!isComplete(ProcessStep.EXTRACT_SOFT_CLIPS) ||
			!isComplete(ProcessStep.EXTRACT_READ_PAIRS)) {
			throw new IllegalStateException("Cannot traverse evidence before evidence extraction");
		}
		if (includeReadPair) {
			final CloseableIterator<SAMRecord> rawPairIt = getContext().getSamReaderIterator(readPair);
			final CloseableIterator<SAMRecord> rawMateIt = getContext().getSamReaderIterator(pairMate);
			final ReadPairEvidenceIterator rawRpIt = new ReadPairEvidenceIterator(this, rawPairIt, rawMateIt);
			getContext().registerBuffer(getSourceFile().getName() + "." + chr + ".rp", rawRpIt);
			final DirectEvidenceWindowedSortingIterator<NonReferenceReadPair> sortedRpIt = new DirectEvidenceWindowedSortingIterator<NonReferenceReadPair>(getContext(), getReadPairSortWindowSize(), applyBlacklistFilter(rawRpIt));
			getContext().registerBuffer(getSourceFile().getName() + "." + chr + ".rp", sortedRpIt);
			final CloseableIterator<NonReferenceReadPair> finalrpIt = new AutoClosingIterator<NonReferenceReadPair>(sortedRpIt,
					Lists.<Closeable>newArrayList(rawRpIt, rawPairIt, rawMateIt));
			itList.add((CloseableIterator<DirectedEvidence>)(Object)finalrpIt);
		}
		if (includeSoftClip) {
			final List<Closeable> scToClose = Lists.newArrayList();
			final CloseableIterator<SAMRecord> rawScIt = getContext().getSamReaderIterator(softClip);
			scToClose.add(rawScIt);
			CloseableIterator<SoftClipEvidence> scIt = new SoftClipEvidenceIterator(this, rawScIt);
			scToClose.add(scIt);
			if (isRealignmentComplete(false)) {
				//log.debug("Realignment is complete for ", input);
				CloseableIterator<SAMRecord> rawRealignIt = getContext().getSamReaderIterator(realigned);
				scToClose.add(rawRealignIt);
				scIt = new RealignedSoftClipEvidenceIterator(scIt, rawRealignIt);
				getContext().registerBuffer(getSourceFile().getName() + "." + chr + ".sc.ra", (RealignedSoftClipEvidenceIterator)scIt);
				scToClose.add(scIt);
			} else {
				//log.info("Realigned soft clip evidence not present due to missing realignment bam ", realigned);
			}
			// sort into evidence order
			DirectEvidenceWindowedSortingIterator<SoftClipEvidence> dit = new DirectEvidenceWindowedSortingIterator<SoftClipEvidence>(getContext(), getSoftClipSortWindowSize(), applyBlacklistFilter(scIt));
			getContext().registerBuffer(getSourceFile().getName() + "." + chr + ".sc", dit);
			scIt =  new AutoClosingIterator<SoftClipEvidence>(dit, scToClose);
			itList.add((CloseableIterator<DirectedEvidence>)(Object)scIt);
		}
		if (includeSoftClipRemote) {
			if (!isRealignmentComplete(false)) {
				throw new IllegalArgumentException(String.format("Realignment resorting not complete for ", input));
			}
			CloseableIterator<SAMRecord> rsrRawIt = getContext().getSamReaderIterator(remoteRealigned);
			CloseableIterator<SAMRecord> sssRawItf = getContext().getSamReaderIterator(remoteSoftClip);
			CloseableIterator<SAMRecord> sssRawItb = getContext().getSamReaderIterator(remoteSoftClip);
			RealignedRemoteSoftClipEvidenceIterator remoteScIt = new RealignedRemoteSoftClipEvidenceIterator(this, rsrRawIt, sssRawItf, sssRawItb);
			getContext().registerBuffer(getSourceFile().getName() + "." + chr + ".scr", remoteScIt);
			DirectEvidenceWindowedSortingIterator<RealignedRemoteSoftClipEvidence> dit = new DirectEvidenceWindowedSortingIterator<RealignedRemoteSoftClipEvidence>(getContext(), getSoftClipSortWindowSize(), applyBlacklistFilter(remoteScIt));
			getContext().registerBuffer(getSourceFile().getName() + "." + chr + ".scr", dit);
			CloseableIterator<RealignedRemoteSoftClipEvidence> sortedRemoteScIt = new AutoClosingIterator<RealignedRemoteSoftClipEvidence>(dit, ImmutableList.<Closeable>of(remoteScIt, rsrRawIt, sssRawItf, sssRawItb));
			itList.add((CloseableIterator<DirectedEvidence>)(Object)sortedRemoteScIt);
		}
		CloseableIterator<DirectedEvidence> mergedIt = new AutoClosingMergedIterator<DirectedEvidence>(itList, DirectedEvidenceOrder.ByNatural);
		if (Defaults.SANITY_CHECK_ITERATORS) {
			mergedIt = new AutoClosingIterator<DirectedEvidence>(new OrderAssertingIterator<DirectedEvidence>(mergedIt, DirectedEvidenceOrder.ByNatural), ImmutableList.<Closeable>of(mergedIt));
		}
		return new AsyncBufferedIterator<DirectedEvidence>(mergedIt, input.getName() + "-" + chr, getContext().getConfig().async_bufferCount, getContext().getConfig().async_bufferSize);
	}
	private <T extends DirectedEvidence> IntervalBedFilteringIterator<T> applyBlacklistFilter(Iterator<T> it) {
		return new IntervalBedFilteringIterator<T>(getBlacklistedRegions(), it, getBlacklistFilterMargin());
	}
	private int getBlacklistFilterMargin() {
		return getMaxConcordantFragmentSize();
	}
	private int getSoftClipSortWindowSize() {
		// worst case: forward with small clip followed by large homologous clip 
		// MMMMMMMMMM>
		//          <<<<<<<<<<M with microhomology the full length of the read
		return 2 * Math.max(getMaxReadMappedLength(), getMaxReadLength()) + 1;
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
					if (getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH != null &&
							getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH >= 100000 &&  
							getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH >= getMetrics().getInsertSizeMetrics().MAX_INSERT_SIZE / 2 &&
							getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH >= 10 * getMaxReadLength()) {
						String msg = String.format("Proper pair flag indicates fragment size of %d is expected!"
								+ " Realign with an aligner that consider fragment size when setting the proper pair flag or "
								+ " specify fixed or percentage bounds for read pair concordance for %s."
								+ " Insert size distribution counts can be found in %s",
								getMetrics().getIdsvMetrics().MAX_PROPER_PAIR_FRAGMENT_LENGTH,
								input,
								getContext().getFileSystemContext().getIdsvMetrics(input)); 
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
	@Override
	public int getMinConcordantFragmentSize() {
		int min = Math.max(getMaxReadLength(), getReadPairConcordanceCalculator().minConcordantFragmentSize());
		return min;
	}
	public int getMaxReadLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_LENGTH;
	}
	public int getMaxReadMappedLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_MAPPED_LENGTH;
	}
	public int getExpectedFragmentSize() {
		return Math.min((int)getMetrics().getInsertSizeMetrics().MEDIAN_INSERT_SIZE, getMaxConcordantFragmentSize());
	}
	public ProcessingContext getProcessContext() {
		return getContext();
	}
	public IntervalBed getBlacklistedRegions() {
		if (blacklist == null) {
			if (!isComplete(ProcessStep.CALCULATE_METRICS)) {
				completeSteps(EnumSet.of(ProcessStep.CALCULATE_METRICS));
			}
			try {
				blacklist = new IntervalBed(
						getContext().getDictionary(),
						getContext().getLinear(),
						getContext().getFileSystemContext().getCoverageBlacklistBed(getSourceFile()));
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			blacklist = IntervalBed.merge(getContext().getDictionary(), getContext().getLinear(), ImmutableList.of(blacklist, getContext().getBlacklistedRegions()));
		}
		return blacklist;
	}
}