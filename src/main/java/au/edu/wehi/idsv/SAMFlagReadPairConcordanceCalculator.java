package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.metrics.IdsvMetrics;

public class SAMFlagReadPairConcordanceCalculator extends ReadPairConcordanceCalculator {
	private final IdsvMetrics metrics;
	public SAMFlagReadPairConcordanceCalculator(IdsvMetrics metrics) {
		this.metrics = metrics;
	}
	@Override
	public int maxConcordantFragmentSize() {
		return metrics.MAX_PROPER_PAIR_FRAGMENT_LENGTH;
	}
	@Override
	public boolean isConcordant(SAMRecord local) {
		return local.getReadPairedFlag() && local.getProperPairFlag();
	}
}
