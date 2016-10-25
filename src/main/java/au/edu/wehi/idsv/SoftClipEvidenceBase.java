package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import java.util.Arrays;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

public abstract class SoftClipEvidenceBase implements DirectedEvidence {
	private final SAMRecord record;
	private final BreakendSummary location;
	/**
	 * Lazily computed evidence identifier
	 */
	private String evidenceID = null;
	public SoftClipEvidenceBase(BreakendDirection direction, SAMRecord record) {
		this.record = record;
		int pos = direction == BreakendDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart();
		this.location = new BreakendSummary(record.getReferenceIndex(), direction, pos, pos);
	}
	public static int getSoftClipLength(BreakendDirection direction, SAMRecord record) {
		return direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(record) : SAMRecordUtil.getStartSoftClipLength(record); 
	}
	/**
	 * Gets the soft clip evidenceID that would be associated with the given record 
	 * @param direction
	 * @param record
	 * @return evidenceID
	 */
	public static String getEvidenceID(BreakendDirection direction, SAMRecord record) {
		return buildSoftClipEvidenceID(direction, record, null).toString();
	}
	private static StringBuilder buildSoftClipEvidenceID(BreakendDirection direction, SAMRecord record, Integer indelOffset) {
		StringBuilder sb = new StringBuilder();
		sb.append(direction.toChar());
		if (indelOffset != null) {
			sb.append(indelOffset.toString());
			sb.append('_');
		}
		sb.append(record.getReadName());
		if (record.getReadPairedFlag()) {
			if (record.getFirstOfPairFlag()) {
				sb.append(SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX);
			} else {
				sb.append(SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX);
			}
		}
		return sb;
	}
	protected Integer indelOffset() {
		return null;
	}
	private StringBuilder buildEvidenceID() {
		return buildSoftClipEvidenceID(location.direction, record, indelOffset());
	}
	@Override
	public String getEvidenceID() {
		if (evidenceID == null) {
			evidenceID = buildEvidenceID().toString();
		}
		return evidenceID;
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
	public byte[] getAnchorSequence() {
		return getAnchor(record.getReadBases());
	}
	@Override
	public byte[] getAnchorQuality() {
		return getAnchor(record.getBaseQualities());
	}
	private byte[] getAnchor(byte[] data) {
		int len = data.length - getSoftClipLength();
		if (location.direction == BreakendDirection.Forward) {
			return Arrays.copyOfRange(data, 0, len);
		} else {
			return Arrays.copyOfRange(data, data.length - len, data.length);
		}
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public int getSoftClipLength() {
		return getSoftClipLength(location.direction, record); 
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
	@Override
	public int getLocalMapq() {
		return record.getMappingQuality();
	}
	@Override
	public int getLocalBaseLength() {
		return record.getReadLength() - SAMRecordUtil.getStartSoftClipLength(record) - SAMRecordUtil.getEndSoftClipLength(record);
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
		return "SC" + getSoftClipLength() + " " + getBreakendSummary().toString() + " " + getEvidenceID();
	}
	/**
	 * Dovetailing reads do not support SVs, they are caused by fragment size
	 * less than read length
	 * 
	 * =======> <=======
	 * 
	 * 
	 * @param expectedOrientation
	 *            read pair orientation
	 * @return true if the soft clip is due to a fragment size smaller than the
	 *         read length
	 */
	public boolean isDovetailing(int dovetailMargin) {
		if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()
				|| record.getReadUnmappedFlag())
			return false;
		return record.getMateReferenceIndex() == record.getReferenceIndex()
				&& Math.abs(record.getAlignmentStart()
						- record.getMateAlignmentStart()) <= dovetailMargin
				// dovetails happen on the 3' end of the read for FR
				&& ((location.direction == BreakendDirection.Forward && !record
						.getReadNegativeStrandFlag()) || (location.direction == BreakendDirection.Backward && record
						.getReadNegativeStrandFlag()));
	}
	@Override
	public boolean isBreakendExact() {
		return true;
	}
}
