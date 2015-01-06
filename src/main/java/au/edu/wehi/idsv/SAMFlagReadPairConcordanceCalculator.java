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
	public boolean isConcordant(SAMRecord read) {
		return read.getReadPairedFlag() && read.getProperPairFlag();
	}
	@Override
	public boolean isConcordant(SAMRecord read1, SAMRecord read2) {
		return isConcordant(read1);
	}
}
