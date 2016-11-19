package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang.NotImplementedException;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Structural variation evidence based on read pairs from a single SAM/BAM.  
 * @author Daniel Cameron
 *
 */
public class SAMEvidenceSource extends EvidenceSource {
	private final SamReaderFactory factory = SamReaderFactory.makeDefault();
	private final File input;
	private final int sourceCategory;
	private final ReadPairConcordanceMethod rpcMethod;
	private final int rpcMinFragmentSize;
	private final int rpcMaxFragmentSize;
	private final double rpcConcordantPercentage;
	private IdsvSamFileMetrics metrics;
	private IntervalBed blacklist;
	private ReadPairConcordanceCalculator rpcc;
	private List<Function<SAMRecord, SAMRecord>> samTransforms;
	private List<SamRecordFilter> samFilters;
	private List<DirectedEvidence> evidenceFilters;
	public SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory) {
		this(processContext, file, sourceCategory, ReadPairConcordanceMethod.SAM_FLAG, 0, 0, 0);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory, int minFragmentSize, int maxFragmentSize) {
		this(processContext, file, sourceCategory, ReadPairConcordanceMethod.FIXED, minFragmentSize, maxFragmentSize, 0);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory, double concordantPercentage) {
		this(processContext, file, sourceCategory, ReadPairConcordanceMethod.PERCENTAGE, 0, 0, concordantPercentage);
	}
	protected SAMEvidenceSource(ProcessingContext processContext, File file, int sourceCategory,
			ReadPairConcordanceMethod rpcMethod, int rpcMinFragmentSize, int rpcMaxFragmentSize, double rpcConcordantPercentage) {
		super(processContext, file);
		this.input = file;
		this.sourceCategory = sourceCategory;
		this.rpcMethod = rpcMethod;
		this.rpcMinFragmentSize = rpcMinFragmentSize;
		this.rpcMaxFragmentSize = rpcMaxFragmentSize;
		this.rpcConcordantPercentage = rpcConcordantPercentage;
	}
	public IdsvSamFileMetrics getMetrics() {
		if (metrics == null) {
			File idsvFile = getContext().getFileSystemContext().getIdsvMetrics(input);
			File cigarFile = getContext().getFileSystemContext().getCigarMetrics(input);
			File isFile = getContext().getFileSystemContext().getInsertSizeMetrics(input);
			File mapqFile = getContext().getFileSystemContext().getMapqMetrics(input);
			if (idsvFile.exists() && cigarFile.exists() &&  isFile.exists() && mapqFile.exists()) {
				metrics = new IdsvSamFileMetrics(getContext(), getFile());
			}
		}
		return metrics;
	}
	/**
	 * Grouping category of source. For tumour/normal analysis, normal is considered group 0, tumour group 1
	 * @return
	 */
	public int getSourceCategory() {
		return sourceCategory;
	}
	public CloseableIterator<DirectedEvidence> iterator(final QueryInterval intervals) {
		SamReader reader = factory.open(input);
		CloseableIterator<SAMRecord> it = reader.queryOverlapping(new QueryInterval[] {intervals});
		Iterator<DirectedEvidence> eit = asEvidence(reader, it);
		// TODO: filter DirectedEvidence to within one off the query intervals
		// - use linear genomic coordinates conversion
		throw new NotImplementedException();
	}
	public CloseableIterator<DirectedEvidence> iterator() {
		SamReader reader = factory.open(input);
		CloseableIterator<SAMRecord> it = reader.iterator();
		return asEvidence(reader, it);
	}
	private CloseableIterator<DirectedEvidence> asEvidence(SamReader reader, CloseableIterator<SAMRecord> rawit) {
		Iterator<SAMRecord> it = new AutoClosingIterator<SAMRecord>(rawit, ImmutableList.of(reader));
		
		// TODO: apply SAMRecord transforms (eg mapq, entropy, adapter)
		// TODO: apply SAMRecord filters (eg blacklist)
		
		Iterator<DirectedEvidence> eit = new DirectedEvidenceIterator(it, this);
		eit = applyBlacklistFilter(eit);
		// TODO: apply evidence filters
		// sc length
		// indel length
		// entropy
		eit = new DirectEvidenceWindowedSortingIterator<DirectedEvidence>(getContext(), getSortWindowSize(), eit);
		if (Defaults.SANITY_CHECK_ITERATORS) {
			eit = new PairedEvidenceTracker<DirectedEvidence>(getFile().getName(), eit);
		}
		throw new NotImplementedException();
	}
	private <T extends DirectedEvidence> IntervalBedFilteringIterator<T> applyBlacklistFilter(Iterator<T> it) {
		return new IntervalBedFilteringIterator<T>(getBlacklistedRegions(), it, getBlacklistFilterMargin());
	}
	private int getBlacklistFilterMargin() {
		return getMaxConcordantFragmentSize();
	}
	private int getSortWindowSize() {
		// Soft clip:
		// worst case: forward with small clip followed by large homologous clip 
		// MMMMMMMMMM>
		//          <<<<<<<<<<M with microhomology the full length of the read
		int softClipSortWindowSize = 2 * Math.max(getMaxReadMappedLength(), getMaxReadLength()) + 1;
		// Read pair:
		// worst case scenario is a fully mapped forward read followed by a large soft-clipped backward read at max distance
		// MMMMMMMMMMMMMMMMMM>
		// ^--read length --^
		// ^--------max concordant fragment size-----^
		//                   |-----breakend call-----|
		//                   |----------breakend call--------|
		//                                                   <SSSSSSSSSM
		//                   ^--------max concordant fragment size-----^
		// ^ alignment start                                           ^ alignment start
		int readPairSortWindowSize = getMaxConcordantFragmentSize() + getMaxReadMappedLength() + 1;
		int windowSize = Math.max(softClipSortWindowSize, readPairSortWindowSize);
		return windowSize;
	}
	public ReadPairConcordanceCalculator getReadPairConcordanceCalculator() {
		if (rpcc == null) {
			rpcc = ReadPairConcordanceCalculator.create(rpcMethod, rpcMinFragmentSize, rpcMaxFragmentSize, rpcConcordantPercentage, getMetrics().getInsertSizeDistribution(), getMetrics().getIdsvMetrics());
		}
		return rpcc;
	}
	@Override
	public int getMaxConcordantFragmentSize() {
		int fs = getMaxReadLength();
		fs = Math.max(getMaxReadMappedLength(), fs);
		if (getReadPairConcordanceCalculator() != null) {
			fs = Math.max(getReadPairConcordanceCalculator().maxConcordantFragmentSize(), fs);
		}
		return fs;
	}
	@Override
	public int getMinConcordantFragmentSize() {
		int fs  = getMaxReadLength();
		if (getReadPairConcordanceCalculator() != null) {
			fs = Math.max(getReadPairConcordanceCalculator().minConcordantFragmentSize(), fs);
		}
		return fs;
	}
	public int getMaxReadLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_LENGTH;
	}
	public int getMaxReadMappedLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_MAPPED_LENGTH;
	}
	public int getExpectedFragmentSize() {
		if (getMetrics().getInsertSizeMetrics() != null) {
			return Math.min((int)getMetrics().getInsertSizeMetrics().MEDIAN_INSERT_SIZE, getMaxConcordantFragmentSize());
		}
		return getMaxConcordantFragmentSize();
	}
	public GenomicProcessingContext getProcessContext() {
		return getContext();
	}
	public IntervalBed getBlacklistedRegions() {
		if (blacklist == null) {
			try {
				blacklist = new IntervalBed(
						getContext().getDictionary(),
						getContext().getLinear(),
						getContext().getFileSystemContext().getCoverageBlacklistBed(getFile()));
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			blacklist = IntervalBed.merge(getContext().getDictionary(), getContext().getLinear(), ImmutableList.of(blacklist, getContext().getBlacklistedRegions()));
		}
		return blacklist;
	}
	public static CloseableIterator<DirectedEvidence> mergedIterator(final List<SAMEvidenceSource> source) {
		List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
		for (SAMEvidenceSource bam : source) {
			CloseableIterator<DirectedEvidence> it = bam.iterator();
			toMerge.add(it);
		}
		CloseableIterator<DirectedEvidence> merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge, DirectedEvidenceOrder.ByNatural);
		return merged;
	}
	public static CloseableIterator<DirectedEvidence> mergedIterator(final List<SAMEvidenceSource> source, final QueryInterval intervals) {
		List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
		for (SAMEvidenceSource bam : source) {
			CloseableIterator<DirectedEvidence> it = bam.iterator(intervals);
			toMerge.add(it);
		}
		CloseableIterator<DirectedEvidence> merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge, DirectedEvidenceOrder.ByNatural);
		return merged;
	}
}