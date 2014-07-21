package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public class RealignedSoftClipEvidence extends SoftClipEvidence implements DirectedBreakpoint {
	private final RealignedBreakpoint realignedBreakpoint;
	private final SAMRecord realigned;
	private BreakpointSummary location;
	protected RealignedSoftClipEvidence(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record, SAMRecord realigned) {
		super(processContext, source, direction, record);
		this.realigned = realigned;
		int pos = direction == BreakendDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart();
		BreakendSummary local = new BreakendSummary(record.getReferenceIndex(), direction, pos, pos);
		this.realignedBreakpoint = new RealignedBreakpoint(processContext, local, record.getReadBases(), realigned);
		this.location = this.realignedBreakpoint.getBreakpointSummary(); 
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return location;
	}
	@Override
	public String getUntemplatedSequence() {
		return realignedBreakpoint.getInsertedSequence();
	}
	@Override
	public int getRemoteMapq() {
		return realigned.getMappingQuality();
	}
	@Override
	public int getRemoteBaseLength() {
		return realigned.getReadLength();
	}
	@Override
	public int getRemoteBaseCount() {
		return getRemoteBaseLength();
	}
	@Override
	public int getRemoteMaxBaseQual() {
		return SAMRecordUtil.getMaxReferenceBaseQual(realigned);
	}
	@Override
	public int getRemoteTotalBaseQual() {
		return SAMRecordUtil.getTotalReferenceBaseQual(realigned);
	}
}
