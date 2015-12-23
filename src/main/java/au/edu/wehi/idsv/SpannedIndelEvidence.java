package au.edu.wehi.idsv;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SpannedIndelEvidence extends RealignedSoftClipEvidence {
	private final int indelOffset;
	public SpannedIndelEvidence(SAMEvidenceSource source, BreakendDirection direction, SAMRecord record, SAMRecord realigned, int indelOffset) {
		super(source, direction, record, realigned);
		this.indelOffset = indelOffset;
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
		return (float)getEvidenceSource().getContext().getConfig().getVariantCalling().getModel().scoreIndel(getEvidenceSource().getMetrics(), op, size, getLocalMapq()) / 2;
	}
	@Override
	public float getBreakpointQual() {
		return getBreakendQual();
	}
}
