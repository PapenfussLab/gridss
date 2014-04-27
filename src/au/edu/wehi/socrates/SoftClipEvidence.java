package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;
import au.edu.wehi.socrates.util.SAMRecordSummary;

public class SoftClipEvidence implements DirectedBreakpoint {
	private final SAMRecord record;
	private final SAMRecordSummary summary;
	private final SAMRecord realigned;
	private final BreakpointLocation location;
	private SoftClipEvidence(BreakpointDirection direction, SAMRecord record, SAMRecordSummary summary, SAMRecord realigned) {
		this.record = record;
		this.summary = summary;
		this.realigned = realigned;
		if (realigned == null) {
			this.location = calculateLocation(direction, record, summary);
		} else {
			this.location = calculateInterval(direction, record, summary, realigned);
		}
	}
	public SoftClipEvidence(BreakpointDirection direction, SAMRecord record) {
		this(direction, record, new SAMRecordSummary(record), null);
	}
	public SoftClipEvidence(BreakpointDirection direction, SAMRecord record, SAMRecord realigned) {
		this(direction, record, new SAMRecordSummary(record), realigned);
	}
	public SoftClipEvidence(SoftClipEvidence evidence, SAMRecord realigned) {
		this(evidence.location.direction, evidence.record, evidence.summary, realigned);
	}
	private static double calculateSoftClipQuality(BreakpointDirection direction, SAMRecord record, SAMRecordSummary summary) {
		// TODO: proper quality metrics!
		if (direction == BreakpointDirection.Forward) {
			return summary.getAvgTailClipQuality();
			//return SAMRecordSummary.getEndSoftClipLength(record);
		} else {
			return summary.getAvgHeadClipQuality();
			//return SAMRecordSummary.getStartSoftClipLength(record);
		}
	}
	private static double calculateSoftClipQuality(BreakpointDirection direction, SAMRecord record, SAMRecordSummary summary, SAMRecord realigned) {
		// TOOD: proper quality metrics!
		return Math.min(record.getMappingQuality(), realigned.getMappingQuality());
	}
	private static BreakpointLocation calculateLocation(BreakpointDirection direction, SAMRecord record, SAMRecordSummary summary) {
		int pos = direction == BreakpointDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart(); 
		return new BreakpointLocation(record.getReferenceIndex(), direction, pos, pos, calculateSoftClipQuality(direction, record, summary));
	}
	private static BreakpointLocation calculateRealignedLocation(BreakpointDirection direction, SAMRecord realigned) {
		int targetPosition;
		BreakpointDirection targetDirection;
		if ((direction == BreakpointDirection.Forward && realigned.getReadNegativeStrandFlag()) ||
				(direction == BreakpointDirection.Backward && !realigned.getReadNegativeStrandFlag())) {
			targetDirection = BreakpointDirection.Forward;
			targetPosition = realigned.getAlignmentEnd();
		} else  {
			targetDirection = BreakpointDirection.Backward;
			targetPosition = realigned.getAlignmentStart();
		}
		return new BreakpointLocation(realigned.getReferenceIndex(), targetDirection, targetPosition, targetPosition, realigned.getMappingQuality());
	}
	private static BreakpointLocation calculateInterval(BreakpointDirection direction, SAMRecord record, SAMRecordSummary summary, SAMRecord realigned) {
		BreakpointLocation location = calculateLocation(direction, record, summary);
		if (realigned.getReadUnmappedFlag()) return location;
		return new BreakpointInterval(location, calculateRealignedLocation(direction, realigned),
				calculateSoftClipQuality(direction, record, summary, realigned));
	}
	@Override
	public String getEvidenceID() {
		// need read name, breakpoint direction & which read in pair
		String readNumber = record.getReadPairedFlag() ? record.getFirstOfPairFlag() ? "/1" : "/2" : "";
		return String.format("%s%s%s", record.getReadName(), readNumber, location.direction == BreakpointDirection.Forward ? "f" : "b");
	}
	@Override
	public byte[] getBreakpointSequence() {
		return location.direction == BreakpointDirection.Forward ? summary.getTailClipSequence() : summary.getHeadClipSequence();
	}
	@Override
	public byte[] getBreakpointQuality() {
		return location.direction == BreakpointDirection.Forward ? summary.getTailClipQuality() : summary.getHeadClipQuality();
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public SAMRecord getSoftClipRealignmentSAMRecord() {
		return this.realigned;
	}
	public int getSoftClipLength() {
		return getBreakpointSequence().length;
	}
	public float getAlignedPercentIdentity() {
		return summary.getAlignedPercentIdentity();
	}
	public float getAverageClipQuality() {
		return location.direction == BreakpointDirection.Forward ? summary.getAvgTailClipQuality() : summary.getAvgHeadClipQuality();
	}
	@Override
	public BreakpointLocation getBreakpointLocation() {
		return location;
	}
}
