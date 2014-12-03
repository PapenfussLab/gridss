package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
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
}
