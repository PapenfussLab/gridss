package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SplitIndel;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SpannedIndelEvidence extends RealignedSoftClipEvidence {
	private final int indelOffset;
	private SpannedIndelEvidence remote;
	public static SpannedIndelEvidence create(SAMEvidenceSource source, SAMRecord record, SplitIndel si, int indelOffset) {
		si = si.onPositive();
		SpannedIndelEvidence left = new SpannedIndelEvidence(source, BreakendDirection.Forward, si.leftAnchored, si.leftRealigned, indelOffset);
		SpannedIndelEvidence right = new SpannedIndelEvidence(source, BreakendDirection.Backward, si.rightAnchored, si.rightRealigned, indelOffset);
		left.remote = right;
		right.remote = left;
		return left;
	}
	private SpannedIndelEvidence(SAMEvidenceSource source, BreakendDirection direction, SAMRecord record, SAMRecord realigned, int indelOffset) {
		super(source, direction, record, realigned);
		this.indelOffset = indelOffset;
		assert(getBreakendSummary().direction == direction);
		assert(getBreakendSummary().direction2 == direction.reverse());
	}
	@Override
	protected Integer indelOffset() {
		return indelOffset;
	}
	@Override
	public float getBreakendQual() {
		CigarOperator op = CigarOperator.INSERTION;
		int size = getUntemplatedSequence().length();
		int deletionSize = Math.abs(getBreakendSummary().start2 - getBreakendSummary().start) - 1;
		if (deletionSize > size) {
			op = CigarOperator.DELETION;
			size = deletionSize;
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreIndel(getEvidenceSource().getMetrics(), op, size, getLocalMapq()) / 2;
	}
	@Override
	public float getBreakpointQual() {
		return getBreakendQual();
	}
	@Override
	public SpannedIndelEvidence asRemote() {
		return remote;
	}
}
