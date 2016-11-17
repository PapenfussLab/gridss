package au.edu.wehi.idsv;

import gridss.analysis.InsertSizeDistribution;

public class PercentageReadPairConcordanceCalculator extends FixedSizeReadPairConcordanceCalculator {
	public PercentageReadPairConcordanceCalculator(InsertSizeDistribution isd, double concordantPercent) {
		super(bounds(isd, concordantPercent)[0], bounds(isd, concordantPercent)[1]);
	}
	private static int[] bounds(InsertSizeDistribution isd, double concordantPercent) {
		if (concordantPercent < 0 || concordantPercent > 1) throw new IllegalArgumentException("concordantPercent must be between 0.0 and 1.0");
		double halfDistributionConcordantPercentage = (1.0 - concordantPercent) / 2;
		int minBound = isd.inverseCumulativeProbability(halfDistributionConcordantPercentage);
		int maxBound = isd.inverseCumulativeProbability(1.0 - halfDistributionConcordantPercentage);
		return new int[] { minBound, maxBound };
	}
}
