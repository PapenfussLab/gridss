package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.MathUtil;
import htsjdk.samtools.SAMRecord;

public class SoftClipEvidence extends SingleReadEvidence {
	private final int clipLength;
	protected SoftClipEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd) {
		super(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd);
		this.clipLength = offsetUnmappedEnd - offsetUnmappedStart;
		if (clipLength <= 0) throw new IllegalArgumentException("Read must be soft clipped");
	}
	public static SoftClipEvidence create(SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE) throw new IllegalArgumentException("Missing read bases");
		if (direction == BreakendDirection.Backward) {
			int clipLength = SAMRecordUtil.getStartSoftClipLength(record);
			int offsetLocalStart = clipLength;
			int offsetLocalEnd = record.getReadLength() - (INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : SAMRecordUtil.getEndSoftClipLength(record));
			int offsetUnmappedStart = 0;
			int offsetUnmappedEnd = clipLength;
			int unanchoredWidth = UnanchoredReadUtil.widthOfImprecision(record.getCigar());
			BreakendSummary bs = new BreakendSummary(
					record.getReferenceIndex(),
					direction,
					MathUtil.average(record.getAlignmentStart(), record.getAlignmentStart() + unanchoredWidth),
					record.getAlignmentStart(),
					record.getAlignmentStart() + unanchoredWidth);
			return new SoftClipEvidence(source, record, bs, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd);
		} else {
			int clipLength = SAMRecordUtil.getEndSoftClipLength(record);
			int unanchoredWidth = UnanchoredReadUtil.widthOfImprecision(record.getCigar());
			return new SoftClipEvidence(source, record,
					new BreakendSummary(record.getReferenceIndex(), direction, MathUtil.average(record.getAlignmentEnd() - unanchoredWidth, record.getAlignmentEnd()), record.getAlignmentEnd() - unanchoredWidth, record.getAlignmentEnd()),
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
		sb.append(getBreakendSummary().direction == BreakendDirection.Forward ? "(f)" : "(b)");
	}
}
