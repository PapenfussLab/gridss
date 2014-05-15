package au.edu.wehi.socrates;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.SequenceUtil;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.socrates.vcf.EvidenceAttributes;

public class SoftClipEvidence implements DirectedBreakpoint {
	private final ProcessingContext processContext;
	private final SAMRecord record;
	private final SAMRecord realigned;
	private final BreakendSummary location;
	public SoftClipEvidence(ProcessingContext processContext, BreakendDirection direction, SAMRecord record, SAMRecord realigned) {
		if (record == null) throw new IllegalArgumentException("record is null");
		if (direction == null) throw new IllegalArgumentException("direction is null");
		if (record.getReadUnmappedFlag()) throw new IllegalArgumentException(String.format("record %s is unmapped", record.getReadName()));
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE ) throw new IllegalArgumentException(String.format("record %s missing sequence information", record.getReadName()));
		this.processContext = processContext;
		this.record = record;
		this.realigned = realigned;
		this.location = calculateLocation(direction);
		addSoftClipMetrics();
		if (getSoftClipLength() == 0) {
			throw new IllegalArgumentException(String.format("record %s is not %s soft clipped", record.getReadName(), direction));
		}
	}
	public SoftClipEvidence(ProcessingContext processContext, BreakendDirection direction, SAMRecord record) {
		this(processContext, direction, record, null);
	}
	public SoftClipEvidence(SoftClipEvidence evidence, SAMRecord realigned) {
		this(evidence.processContext, evidence.location.direction, evidence.record, realigned);
	}
	private void addSoftClipMetrics() {
		location.evidence.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 1);
		location.evidence.set(EvidenceAttributes.SOFT_CLIP_MAX_LENGTH, getSoftClipLength());
		location.evidence.set(EvidenceAttributes.SOFT_CLIP_TOTAL_LENGTH, getSoftClipLength());
		if (realigned != null && !realigned.getReadUnmappedFlag()) {
			location.evidence.set(EvidenceAttributes.REALIGN_MAX_LENGTH, getSoftClipLength() - getUntemplatedSequenceLength());
			location.evidence.set(EvidenceAttributes.REALIGN_TOTAL_LENGTH, getSoftClipLength() - getUntemplatedSequenceLength());
			location.evidence.set(EvidenceAttributes.REALIGN_MAX_MAPQ, realigned.getMappingQuality());
			location.evidence.set(EvidenceAttributes.REALIGN_TOTAL_MAPQ, realigned.getMappingQuality());
		}
	}
	private BreakendSummary calculateLocation(BreakendDirection direction) {
		int pos = direction == BreakendDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart(); 
		BreakendSummary location = new BreakendSummary(record.getReferenceIndex(), direction, pos, pos, new EvidenceMetrics());
		if (realigned != null && !realigned.getReadUnmappedFlag()) {
			int targetPosition;
			BreakendDirection targetDirection;
			if ((direction == BreakendDirection.Forward && realigned.getReadNegativeStrandFlag()) ||
					(direction == BreakendDirection.Backward && !realigned.getReadNegativeStrandFlag())) {
				targetDirection = BreakendDirection.Forward;
				targetPosition = realigned.getAlignmentEnd();
			} else  {
				targetDirection = BreakendDirection.Backward;
				targetPosition = realigned.getAlignmentStart();
			}
			BreakendSummary targetLocation = new BreakendSummary(realigned.getReferenceIndex(), targetDirection, targetPosition, targetPosition, null);
			location = new BreakpointSummary(location, targetLocation, new EvidenceMetrics());
		}
		return location;
	}
	@Override
	public String getEvidenceID() {
		// need read name, breakpoint direction & which read in pair
		String readNumber = record.getReadPairedFlag() ? record.getFirstOfPairFlag() ? "/1" : "/2" : "";
		return String.format("%s%s%s", location.direction == BreakendDirection.Forward ? "f" : "b", record.getReadName(), readNumber);
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		return location;
	}
	@Override
	public byte[] getBreakpointSequence() {
		byte[] seq = record.getReadBases();
		if (seq == null) return null;
		if (location.direction == BreakendDirection.Forward) {
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
		if (location.direction == BreakendDirection.Forward) {
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
		return getSoftClipLength(location.direction, record); 
	}
	public static int getSoftClipLength(BreakendDirection direction, SAMRecord record) {
		return direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(record) : SAMRecordUtil.getStartSoftClipLength(record); 
	}
	/**
	 * Number of unmapped bases at the breakpoint 
	 * @return Number of unmapped bases at the breakpoint
	 */
	public int getUntemplatedSequenceLength() {
		if (location instanceof BreakpointSummary) {
			BreakpointSummary interval = (BreakpointSummary)location;
			if (interval.direction2 == BreakendDirection.Forward) {
				return SAMRecordUtil.getEndSoftClipLength(realigned);
			} else {
				return SAMRecordUtil.getStartSoftClipLength(realigned);
			}
		} else {
			return getSoftClipLength();
		}
	}
	/**
	 * Sequence of untemplated bases
	 * @return
	 */
	public String getUntemplatedSequence() {
		int untemplatedSequenceLength = getUntemplatedSequenceLength();
		String s = new String(getBreakpointSequence(), StandardCharsets.US_ASCII);
		if (getBreakendSummary().direction == BreakendDirection.Forward) {
			return s.substring(0, untemplatedSequenceLength);
		} else {
			return s.substring(s.length() - untemplatedSequenceLength);
		}
	}
	/**
	 * 0-100 scaled percentage identity of mapped read bases.
	 * @return percentage of reference-aligned bases that match reference 
	 */
	public float getAlignedPercentIdentity() {
		// final byte[] referenceBases = refSeq.get(sequenceDictionary.getSequenceIndex(rec.getReferenceName())).getBases();
        // rec.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(rec, referenceBases, 0, bisulfiteSequence));
        //if (rec.getBaseQualities() != SAMRecord.NULL_QUALS) {
        // rec.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(rec, referenceBases, 0, bisulfiteSequence));
        Integer nm = record.getIntegerAttribute(SAMTag.NM.name());
		if (nm != null) {
			int refBasesToConsider = record.getReadLength() - SAMRecordUtil.getStartSoftClipLength(record) - SAMRecordUtil.getEndSoftClipLength(record); 
			int refBaseMatches = refBasesToConsider - nm + SequenceUtil.countInsertedBases(record) + SequenceUtil.countDeletedBases(record); 
			return 100.0f * refBaseMatches / (float)refBasesToConsider;
		}
		String md = record.getStringAttribute(SAMTag.MD.name());
		if (StringUtils.isNotEmpty(md)) {
			// Socrates handles this: should we? Which aligners write MD but not NM?
			throw new RuntimeException("Sanity Check Failure: Not Yet Implemented: calculation from reads with MD tag but not NM tag as per Socrates implementation");
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
