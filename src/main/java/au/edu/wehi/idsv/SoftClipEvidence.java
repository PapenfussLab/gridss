package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

public class SoftClipEvidence extends SingleReadEvidence {
	private final int clipLength;
	private static final Log log = Log.getInstance(SoftClipEvidence.class);
	protected SoftClipEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd,
			int unanchoredWidth,
		   	boolean isInAssemblyAnchor) {
		super(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, unanchoredWidth, isInAssemblyAnchor);
		this.clipLength = offsetUnmappedEnd - offsetUnmappedStart;
		if (clipLength <= 0) throw new IllegalArgumentException("Read must be soft clipped");
	}
	public static SoftClipEvidence create(SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE) throw new IllegalArgumentException("Missing read bases");
		int unanchoredWidth = CigarUtil.widthOfImprecision(record.getCigar());
		SoftClipEvidence sce = null;
		if (direction == BreakendDirection.Backward) {
			int clipLength = SAMRecordUtil.getStartSoftClipLength(record);
			int offsetLocalStart = clipLength;
			int offsetLocalEnd = record.getReadLength() - (INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : SAMRecordUtil.getEndSoftClipLength(record));
			int offsetUnmappedStart = 0;
			int offsetUnmappedEnd = clipLength;
			BreakendSummary bs = new BreakendSummary(record.getReferenceIndex(), direction, record.getAlignmentStart());
			boolean isInAssemblyAnchor = isEntirelyContainedInAssemblyAnchor(record, new ChimericAlignment(record),
					new ChimericAlignment(record.getReferenceName(), record.getAlignmentStart(), record.getReadNegativeStrandFlag(),
							new Cigar(ImmutableList.of(
									new CigarElement(clipLength, CigarOperator.M),
									new CigarElement(record.getReadLength() - clipLength, CigarOperator.S))),
							record.getMappingQuality()));
			sce = new SoftClipEvidence(source, record, bs, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, unanchoredWidth, isInAssemblyAnchor);

		} else {
			int clipLength = SAMRecordUtil.getEndSoftClipLength(record);
			boolean isInAssemblyAnchor = isEntirelyContainedInAssemblyAnchor(record, new ChimericAlignment(record),
					new ChimericAlignment(record.getReferenceName(), record.getAlignmentStart(), record.getReadNegativeStrandFlag(),
							new Cigar(ImmutableList.of(
									new CigarElement(record.getReadLength() - clipLength, CigarOperator.S),
									new CigarElement(clipLength, CigarOperator.M))),
							record.getMappingQuality()));
			sce = new SoftClipEvidence(source, record,
					new BreakendSummary(record.getReferenceIndex(), direction, record.getAlignmentEnd()),
					INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : SAMRecordUtil.getStartSoftClipLength(record), record.getReadLength() - clipLength,
					record.getReadLength() - clipLength, record.getReadLength(),
					unanchoredWidth, isInAssemblyAnchor);
		}
		return sce;
	}
	@Override
	public float getBreakendQual() {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			return scoreAssembly();
		}
		return (float)source.getContext().getConfig().getScoring().getModel().scoreSoftClip(source.getMetrics(), this, clipLength, getLocalMapq());
	}
	private float scoreAssembly() {
		AssemblyAttributes attr = new AssemblyAttributes(getSAMRecord());
		int pos = getBreakendAssemblyContigOffset();
		int rp = attr.getSupportingReadCount(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null, source.getContext());
		double rpq = attr.getSupportingQualScore(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null, source.getContext());
		int sc = attr.getSupportingReadCount(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null, source.getContext());
		double scq = attr.getSupportingQualScore(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null, source.getContext());
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreBreakendAssembly(this,
				rp, rpq,
				sc, scq,
				getLocalMapq());
	}
	@Override
	protected String getUncachedEvidenceID() {
		return source.getContext().getEvidenceIDGenerator().getEvidenceID(this);
	}
	@Override
	public boolean isReference() {
		return false;
	}
}
