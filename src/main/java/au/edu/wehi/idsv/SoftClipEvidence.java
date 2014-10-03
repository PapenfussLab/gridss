package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class SoftClipEvidence implements DirectedEvidence {
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	private final SAMRecord record;
	private final BreakendSummary location;
	public static SoftClipEvidence create(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		return create(processContext, source, direction, record, null);
	}
	public static SoftClipEvidence create(SoftClipEvidence evidence, SAMRecord realigned) {
		return create(evidence.processContext, evidence.source, evidence.location.direction, evidence.record, realigned);
	}
	public static SoftClipEvidence create(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record, SAMRecord realigned) {
		if (record == null) throw new IllegalArgumentException("record is null");
		if (direction == null) throw new IllegalArgumentException("direction is null");
		if (record.getReadUnmappedFlag()) throw new IllegalArgumentException(String.format("record %s is unmapped", record.getReadName()));
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE ) throw new IllegalArgumentException(String.format("record %s missing sequence information", record.getReadName()));
		SoftClipEvidence result = null;
		if (realigned != null && !realigned.getReadUnmappedFlag()) {
			try {
				result = new RealignedSoftClipEvidence(processContext, source, direction, record, realigned);
			} catch (CloneNotSupportedException e) {
				throw new RuntimeException(e);
			}
		} else {
			result = new SoftClipEvidence(processContext, source, direction, record);
		}
		if (result.getSoftClipLength() == 0) {
			throw new IllegalArgumentException(String.format("record %s is not %s soft clipped", record.getReadName(), direction));
		}
		return result;
	}
	protected SoftClipEvidence(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		this.processContext = processContext;
		this.source = source;
		this.record = record;
		int pos = direction == BreakendDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart();
		this.location = new BreakendSummary(record.getReferenceIndex(), direction, pos, pos);
	}
	public static int getSoftClipLength(BreakendDirection direction, SAMRecord record) {
		return direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(record) : SAMRecordUtil.getStartSoftClipLength(record); 
	}
	public static String getEvidenceID(BreakendDirection direction, SAMRecord softClip) {
		// need read name, breakpoint direction & which read in pair
		String readNumber = softClip.getReadPairedFlag() ? softClip.getFirstOfPairFlag() ? "/1" : "/2" : "";
		return String.format("%s%s%s", direction == BreakendDirection.Forward ? "f" : "b", softClip.getReadName(), readNumber);
	}
	@Override
	public String getEvidenceID() {
		return getEvidenceID(location.direction, record);
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		return location;
	}
	@Override
	public byte[] getBreakendSequence() {
		return location.direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBases(record) : SAMRecordUtil.getStartSoftClipBases(record);
	}
	@Override
	public byte[] getBreakendQuality() {
		return location.direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBaseQualities(record) : SAMRecordUtil.getStartSoftClipBaseQualities(record);
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public int getSoftClipLength() {
		return getSoftClipLength(location.direction, record); 
	}
	/**
	 * Sequence of untemplated bases
	 * @return
	 */
	public String getUntemplatedSequence() {
		return new String(getBreakendSequence(), StandardCharsets.US_ASCII);
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
		byte[] qual = getBreakendQuality();
		if (qual == null) return 0;
		for (int i = 0; i < qual.length; i++) {
			total += qual[i]; 
		}
		return total / qual.length;
	}
	public int getMappingQuality() {
		return record.getMappingQuality();
	}
	public SAMEvidenceSource getEvidenceSource() {
		return source;
	}
	@Override
	public int getLocalMapq() {
		return record.getMappingQuality();
	}
	@Override
	public int getLocalBaseLength() {
		return record.getReadLength() - SAMRecordUtil.getStartSoftClipLength(record) - SAMRecordUtil.getEndSoftClipLength(record);
	}
	@Override
	public int getLocalBaseCount() {
		return getLocalBaseLength();
	}
	@Override
	public int getLocalMaxBaseQual() {
		return SAMRecordUtil.getMaxReferenceBaseQual(record);
	}
	@Override
	public int getLocalTotalBaseQual() {
		return SAMRecordUtil.getTotalReferenceBaseQual(record);
	}
	@Override
	public String toString() {
		return "SoftClip len=" + getSoftClipLength() + " " + getBreakendSummary().toString() + " " + getSAMRecord().getReadName();
	}
}
