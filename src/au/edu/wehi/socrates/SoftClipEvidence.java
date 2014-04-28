package au.edu.wehi.socrates;

import java.util.Arrays;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMTag;
import net.sf.samtools.util.SequenceUtil;

import org.apache.commons.lang3.StringUtils;

public class SoftClipEvidence implements DirectedBreakpoint {
	private final SAMRecord record;
	private final SAMRecord realigned;
	private /*final*/ BreakpointLocation location;
	public SoftClipEvidence(BreakpointDirection direction, SAMRecord record, SAMRecord realigned) {
		if (record == null) throw new IllegalArgumentException("record is null");
		if (direction == null) throw new IllegalArgumentException("direction is null");
		this.record = record;
		this.realigned = realigned;
		this.location = calculateLocation(direction);
		// need to defer quality calculation to here since it depends on location information 
		setQuality(calculateSoftClipQuality());
		if (getSoftClipLength() == 0) throw new IllegalArgumentException(String.format("record %s is not %s soft clipped", record.getReadName(), direction));
	}
	public SoftClipEvidence(BreakpointDirection direction, SAMRecord record) {
		this(direction, record, null);
	}
	public SoftClipEvidence(SoftClipEvidence evidence, SAMRecord realigned) {
		this(evidence.location.direction, evidence.record, realigned);
	}
	private double calculateSoftClipQuality() {
		// TOOD: proper quality metrics!
		if (realigned == null) {
			return getAverageClipQuality();
		} else if (realigned.getReadUnmappedFlag()) {
			return getAverageClipQuality();
		} else {
			return Math.min(record.getMappingQuality(), realigned.getMappingQuality());
		}
	}
	private void setQuality(double qual) {
		if (location instanceof BreakpointInterval) {
			location = new BreakpointInterval((BreakpointInterval)location, qual); 
		} else {
			location = new BreakpointLocation(location, qual);
		}
	}
	private BreakpointLocation calculateLocation(BreakpointDirection direction) {
		int pos = direction == BreakpointDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart(); 
		BreakpointLocation location = new BreakpointLocation(record.getReferenceIndex(), direction, pos, pos, 0);
		if (realigned != null && !realigned.getReadUnmappedFlag()) {
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
			BreakpointLocation targetLocation = new BreakpointLocation(realigned.getReferenceIndex(), targetDirection, targetPosition, targetPosition, 0);
			location = new BreakpointInterval(location, targetLocation, 0);
		}
		return location;
	}
	@Override
	public String getEvidenceID() {
		// need read name, breakpoint direction & which read in pair
		String readNumber = record.getReadPairedFlag() ? record.getFirstOfPairFlag() ? "/1" : "/2" : "";
		return String.format("%s%s%s", location.direction == BreakpointDirection.Forward ? "f" : "b", record.getReadName(), readNumber);
	}
	@Override
	public BreakpointLocation getBreakpointLocation() {
		return location;
	}
	@Override
	public byte[] getBreakpointSequence() {
		byte[] seq = record.getReadBases();
		if (seq == null) return null;
		if (location.direction == BreakpointDirection.Forward) {
			seq = Arrays.copyOfRange(seq, record.getReadLength() - getSoftClipLength(), record.getReadLength()); 
		} else {
			seq = Arrays.copyOfRange(seq, 0, getSoftClipLength());
		}
		return seq;
	}
	@Override
	public byte[] getBreakpointQuality() {
		byte[] seq = record.getBaseQualities();
		if (seq == null) return null;
		if (seq == SAMRecord.NULL_QUALS) return null;
		if (location.direction == BreakpointDirection.Forward) {
			seq = Arrays.copyOfRange(seq, record.getReadLength() - getSoftClipLength(), record.getReadLength()); 
		} else {
			seq = Arrays.copyOfRange(seq, 0, getSoftClipLength());
		}
		return seq;
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public SAMRecord getSoftClipRealignmentSAMRecord() {
		return this.realigned;
	}
	public int getSoftClipLength() {
		return this.location.direction == BreakpointDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(record) : SAMRecordUtil.getStartSoftClipLength(record); 
	}
	/**
	 * Number of unmapped bases at the breakpoint 
	 * @return Number of unmapped bases at the breakpoint
	 */
	public int getUntemplatedSequenceLength() {
		if (location instanceof BreakpointInterval) {
			BreakpointInterval interval = (BreakpointInterval)location;
			if (interval.direction2 == BreakpointDirection.Forward) {
				return SAMRecordUtil.getEndSoftClipLength(realigned);
			} else {
				return SAMRecordUtil.getStartSoftClipLength(realigned);
			}
		} else {
			return getSoftClipLength();
		}
	}
	public float getAlignedPercentIdentity() {
		// final byte[] referenceBases = refSeq.get(sequenceDictionary.getSequenceIndex(rec.getReferenceName())).getBases();
        // rec.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(rec, referenceBases, 0, bisulfiteSequence));
        //if (rec.getBaseQualities() != SAMRecord.NULL_QUALS) {
        // rec.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(rec, referenceBases, 0, bisulfiteSequence));
		
        Integer nm = record.getIntegerAttribute(SAMTag.NM.name());
		if (nm != null) {
			return nm - SequenceUtil.countInsertedBases(record) - SequenceUtil.countDeletedBases(record);
		}
		String md = record.getStringAttribute(SAMTag.MD.name());
		if (StringUtils.isNotEmpty(md)) {
			// Socrates handles this: should we? Which aligners write MD but not NM?
			throw new RuntimeException("TODO: Not Yet Implemented: calculate from reads with MD tag but not NM tag as per Socrates implementation. User Workaround: create NM tags for reads.");
		}
		throw new IllegalStateException(String.format("Read %s missing NM tag", record.getReadName()));
	}
	public float getAverageClipQuality() {
		float total = 0;
		byte[] qual = getBreakpointQuality();
		if (qual == null) return 0;
		for (int i = 0; i < qual.length; i++) {
			total += qual[i]; 
		}
		return total / qual.length;
	}
	public int getMappingQuality() {
		return record.getMappingQuality();
	}
}
