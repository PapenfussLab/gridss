package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

import au.edu.wehi.socrates.util.SAMRecordSummary;

public class SoftClipEvidence implements DirectedBreakpoint {
	public static final char BREAKPOINT_ID_CHAR = 'r';
	private BreakpointDirection direction;
	private SAMRecord record;
	private SAMRecordSummary summary;
	private SAMRecord realigned;
	public SoftClipEvidence(BreakpointDirection direction, SAMRecord record) {
		this.direction = direction;
		this.record = record;
		this.summary = new SAMRecordSummary(record);
	}
	@Override
	public String getEvidenceID() {
		return String.format("%c#%d#%s", BREAKPOINT_ID_CHAR, record.getFirstOfPairFlag(), record.getReadName());
	}
	@Override
	public BreakpointDirection getBreakpointDirection() {
		return direction;
	}
	@Override
	public byte[] getBreakpointSequence() {
		return direction == BreakpointDirection.Forward ? summary.getTailClipSequence() : summary.getHeadClipSequence();
	}
	@Override
	public byte[] getBreakpointQuality() {
		return direction == BreakpointDirection.Forward ? summary.getTailClipQuality() : summary.getHeadClipQuality();
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public int getSoftClipLength() {
		return getBreakpointSequence().length;
	}
	public float getAlignedPercentIdentity() {
		return summary.getAlignedPercentIdentity();
	}
	public float getAverageClipQuality() {
		return direction == BreakpointDirection.Forward ? summary.getAvgTailClipQuality() : summary.getAvgHeadClipQuality();
	}
	@Override
	public int getReferenceIndex() {
		return record.getReferenceIndex();
	}
	@Override
	public long getWindowStart() {
		return direction == BreakpointDirection.Forward ? record.getAlignmentEnd() : (record.getAlignmentStart() - 1);
	}
	@Override
	public long getWindowEnd() {
		return getWindowStart();
	}
	public int getMappingQuality() {
		return record.getMappingQuality();
	}
	@Override
	public void setRealigned(SAMRecord realigned) {
		this.realigned = realigned;
	}
	@Override
	public SAMRecord getRealigned() {
		return realigned;
	}
}
