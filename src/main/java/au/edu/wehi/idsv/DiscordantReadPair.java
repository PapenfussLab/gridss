package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class DiscordantReadPair extends NonReferenceReadPair implements DirectedBreakpoint {
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
	public String toString() {
		return String.format("DP %s MQ=%d,%d RN=%s", getBreakendSummary(), getLocalMapq(), getRemoteMapq(), getEvidenceID());
	}
	@Override
	public String getUntemplatedSequence() {
		return "";
	}
	@Override
	public float getBreakendQual() {
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreReadPair(
				getEvidenceSource().getMetrics(),
				SAMRecordUtil.calculateFragmentSize(getLocalledMappedRead(), getNonReferenceRead(), PairOrientation.FR),
				getLocalMapq(),
				Integer.MAX_VALUE);
	}
	@Override
	public float getBreakpointQual() {
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreReadPair(
				getEvidenceSource().getMetrics(),
				SAMRecordUtil.calculateFragmentSize(getLocalledMappedRead(), getNonReferenceRead(), PairOrientation.FR),
				getLocalMapq(),
				getRemoteMapq());
	}
	public DiscordantReadPair asRemote() {
		return (DiscordantReadPair)NonReferenceReadPair.create(getNonReferenceRead(), getLocalledMappedRead(), getEvidenceSource());
	}
	@Override
	public String getHomologySequence() {
		return null;
	}
	@Override
	public int getHomologyAnchoredBaseCount() {
		return 0;
	}
}
