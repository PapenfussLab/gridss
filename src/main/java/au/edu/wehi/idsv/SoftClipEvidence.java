package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

public class SoftClipEvidence extends SingleReadEvidence {
	private final int clipLength;
	protected SoftClipEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd,
			int unanchoredWidth) {
		super(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, unanchoredWidth);
		this.clipLength = offsetUnmappedEnd - offsetUnmappedStart;
		if (clipLength <= 0) throw new IllegalArgumentException("Read must be soft clipped");
	}
	public static SoftClipEvidence create(SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE) throw new IllegalArgumentException("Missing read bases");
		int unanchoredWidth = CigarUtil.widthOfImprecision(record.getCigar());
		SoftClipEvidence sce;
		if (direction == BreakendDirection.Backward) {
			int clipLength = SAMRecordUtil.getStartSoftClipLength(record);
			int offsetLocalStart = clipLength;
			int offsetLocalEnd = record.getReadLength() - (INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : SAMRecordUtil.getEndSoftClipLength(record));
			int offsetUnmappedStart = 0;
			int offsetUnmappedEnd = clipLength;
			BreakendSummary bs = new BreakendSummary(record.getReferenceIndex(), direction, record.getAlignmentStart());
			sce = new SoftClipEvidence(source, record, bs, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, unanchoredWidth);
		} else {
			int clipLength = SAMRecordUtil.getEndSoftClipLength(record);
			sce = new SoftClipEvidence(source, record,
					new BreakendSummary(record.getReferenceIndex(), direction, record.getAlignmentEnd()),
					INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : SAMRecordUtil.getStartSoftClipLength(record), record.getReadLength() - clipLength,
					record.getReadLength() - clipLength, record.getReadLength(),
					unanchoredWidth);
		}
		return sce;
	}
	@Override
	public float getBreakendQual() {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			return scoreAssembly();
		}
		return (float)source.getContext().getConfig().getScoring().getModel().scoreSoftClip(source.getMetrics(), clipLength, getLocalMapq());
	}
	private float scoreAssembly() {
		if (getBreakendSequence().length == 0) return 0;
		AssemblyAttributes attr = new AssemblyAttributes(getSAMRecord());
		int rp = attr.getAssemblySupportCountReadPair();
		double rpq = attr.getAssemblySupportReadPairQualityScore();
		int sc = attr.getAssemblySupportCountSoftClip();
		double scq =  attr.getAssemblySupportSoftClipQualityScore();
		if (source.getContext().getAssemblyParameters().excludeNonSupportingEvidence) {
			rp -= attr.getAssemblyNonSupportingReadPairCount();
			rpq -= attr.getAssemblyNonSupportingReadPairQualityScore();
			sc -= attr.getAssemblyNonSupportingSoftClipCount();
			scq -= attr.getAssemblyNonSupportingSoftClipQualityScore();
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreBreakendAssembly(
				rp, rpq,
				sc, scq,
				getLocalMapq());
	}
	@Override
	protected void buildEvidenceID(StringBuilder sb) {
		super.buildEvidenceID(sb);
		sb.append(getBreakendSummary().direction == BreakendDirection.Forward ? "(f)" : "(b)");
	}
}
