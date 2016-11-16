package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

public class SoftClipReadEvidence extends SingleReadEvidence {
	private final int clipLength;
	protected SoftClipReadEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd) {
		super(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd);
		this.clipLength = offsetUnmappedEnd - offsetUnmappedStart;
		if (clipLength <= 0) throw new IllegalArgumentException("Read must be soft clipped");
	}
	public static SoftClipReadEvidence create(SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE) throw new IllegalArgumentException("Missing read bases");
		if (direction == BreakendDirection.Backward) {
			int clipLength = SAMRecordUtil.getStartSoftClipLength(record);
			return new SoftClipReadEvidence(source, record,
					new BreakendSummary(record.getReferenceIndex(), direction, record.getAlignmentStart()),
					clipLength, record.getReadLength() - (INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : SAMRecordUtil.getEndSoftClipLength(record)),
					0, clipLength);
		} else {
			int clipLength = SAMRecordUtil.getEndSoftClipLength(record);
			return new SoftClipReadEvidence(source, record,
					new BreakendSummary(record.getReferenceIndex(), direction, record.getAlignmentEnd()),
					INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : SAMRecordUtil.getStartSoftClipLength(record), record.getReadLength() - clipLength,
					record.getReadLength() - clipLength, record.getReadLength());
		}
	}
	@Override
	public float getBreakendQual() {
		return (float)source.getContext().getConfig().getScoring().getModel().scoreSoftClip(source.getMetrics(), clipLength, getLocalMapq());
	}
	@Override
	protected void buildEvidenceID(StringBuilder sb) {
		super.buildEvidenceID(sb);
		sb.append(getBreakendSummary().direction == BreakendDirection.Forward ? 'f' : 'b');
	}
}
