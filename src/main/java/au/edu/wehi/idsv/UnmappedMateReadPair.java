package au.edu.wehi.idsv;

import au.edu.wehi.idsv.metrics.IdsvMetrics;
import htsjdk.samtools.SAMRecord;

public class UnmappedMateReadPair extends NonReferenceReadPair {
	protected UnmappedMateReadPair(SAMRecord local, SAMRecord remote, SAMEvidenceSource source, ProcessingContext processContext) {
		super(local, remote, source, processContext);
	}
	@Override
	public String toString() {
		return String.format("UM %s MQ=%d RN=%s", getBreakendSummary(), getLocalMapq(), getLocalledMappedRead().getReadName());
	}
	@Override
	public float getBreakendQual() {
		IdsvMetrics metrics = getEvidenceSource().getMetrics().getIdsvMetrics();
		// completely unmapped read pairs are excluded for consistency with sc and dp calculation
		double score = -10 * Math.log10((double)metrics.READ_PAIRS_ONE_MAPPED / (double)(metrics.READ_PAIRS - metrics.READ_PAIRS_ZERO_MAPPED));
		score = Math.min(score, getLocalMapq());
		return (float)score;
	}
}
