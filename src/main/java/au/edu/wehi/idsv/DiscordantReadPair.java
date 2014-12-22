package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.InsertSizeDistribution;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class DiscordantReadPair extends NonReferenceReadPair implements DirectedBreakpoint {
	protected DiscordantReadPair(SAMRecord local, SAMRecord remote, SAMEvidenceSource source, ProcessingContext processContext) {
		super(local, remote, source, processContext);
		assert(!remote.getReadUnmappedFlag());
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return (BreakpointSummary)super.getBreakendSummary();
	}
	@Override
	public int getRemoteMapq() {
		return getNonReferenceRead().getMappingQuality();
	}
	@Override
	public int getRemoteBaseLength() {
		return getNonReferenceRead().getReadLength();
	}

	@Override
	public int getRemoteBaseCount() {
		return getNonReferenceRead().getReadLength();
	}

	@Override
	public int getRemoteMaxBaseQual() {
		return SAMRecordUtil.getMaxReferenceBaseQual(getNonReferenceRead());
	}

	@Override
	public int getRemoteTotalBaseQual() {
		return SAMRecordUtil.getTotalReferenceBaseQual(getNonReferenceRead());
	}
	@Override
	public String toString() {
		return String.format("DP %s MQ=%d,%d RN=%s", getBreakendSummary(), getLocalMapq(), getRemoteMapq(), getLocalledMappedRead().getReadName());
	}
	@Override
	public String getUntemplatedSequence() {
		return "";
	}
	/**
	 * 
	 * @param source
	 * @param fragmentSize size of fragment, 0 indicates read pair does not support a reference-allele fragment at all
	 * @param localMapq
	 * @param remoteMapq
	 * @return
	 */
	protected static float dpPhred(SAMEvidenceSource source, int fragmentSize, int localMapq, int remoteMapq) {
		double pairsFromFragmentDistribution = 0;
		InsertSizeDistribution isd = source.getMetrics().getInsertSizeDistribution();
		IdsvMetrics metrics = source.getMetrics().getIdsvMetrics();
		if (fragmentSize > 0) {
			if (fragmentSize >= isd.getSupportLowerBound() && fragmentSize <= isd.getSupportUpperBound()) {
				double pr = source.getMetrics().getInsertSizeDistribution().cumulativeProbability(fragmentSize);
				if (pr < 0.5) {
					pr = 1.0 - pr;
				}
				pairsFromFragmentDistribution = pr * isd.getTotalMappedPairs();
			}
		}
		double totalPairs = metrics.READ_PAIRS_BOTH_MAPPED;
		double dpPairs = totalPairs - isd.getTotalMappedPairs() + pairsFromFragmentDistribution;
		double score = -10 * Math.log10(dpPairs / totalPairs);
		score = Math.min(score, localMapq);
		score = Math.min(score, remoteMapq);
		return (float)score;
	}
	@Override
	public float getBreakendQual() {
		return dpPhred(getEvidenceSource(),
				SAMRecordUtil.calculateFragmentSize(getLocalledMappedRead(), getNonReferenceRead()),
				getLocalMapq(),
				Integer.MAX_VALUE);
	}
	@Override
	public float getBreakpointQual() {
		return dpPhred(getEvidenceSource(),
				SAMRecordUtil.calculateFragmentSize(getLocalledMappedRead(), getNonReferenceRead()),
				getLocalMapq(),
				getRemoteMapq());
	}
}
