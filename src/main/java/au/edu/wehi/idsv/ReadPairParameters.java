package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

public class ReadPairParameters {
	/**
	 * Use SAM flag to extract read pairs
	 */
	public boolean useProperPairFlag = true;
	/**
	 * Percentage (0.0-1.0) of reads considered concordant and ignored for read pair analysis.
	 */
	public double concordantPercent = 0;
	/**
	 * Minimum MAPQ of local read anchor to considered as evidence
	 */
	public int minLocalMapq = 5;
	private double getCordantPercentageBound() {
		double halfDistributionConcordantPercentage = concordantPercent + (1.0 - concordantPercent) / 2; // 0.95 + (1 - 0.95) / 2 = 0.975
		return halfDistributionConcordantPercentage;
	}
	public double getCordantPercentageUpperBound() {
		return getCordantPercentageBound();
	}
	public double getCordantPercentageLowerBound() {
		return 1.0 - getCordantPercentageBound();
	}
	public boolean meetsEvidenceCritera(NonReferenceReadPair evidence) {
		return meetsEvidenceCritera(evidence.getLocalledMappedRead());
	}
	public boolean meetsEvidenceCritera(SAMRecord local) {
		return local.getMappingQuality() >= minLocalMapq;
	}
}
