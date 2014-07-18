package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public class DiscordantReadPair extends NonReferenceReadPair implements DirectedBreakpoint {
	@Override
	public float getPhredLogLikelihoodRatio() {
		// TODO: real metrics
		return Math.min(getLocalMapq(), getRemoteMapq());
	}
	protected DiscordantReadPair(SAMRecord local, SAMRecord remote, SAMEvidenceSource source) {
		super(local, remote, source);
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
}
