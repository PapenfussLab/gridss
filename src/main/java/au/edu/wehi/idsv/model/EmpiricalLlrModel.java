package au.edu.wehi.idsv.model;

import htsjdk.samtools.CigarOperator;
import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.metrics.InsertSizeDistribution;

/**
 * Variant scoring model that scores variants according to the library empirical
 * distribution of the relevant feature
 *  
 * @author Daniel Cameron
 *
 */
public class EmpiricalLlrModel implements VariantScoringModel {
	
	private double llr(double prEgivenMR, double prEgivenMV, double prM) {
		double likelihoodRatio =  (prEgivenMR + prM * (prEgivenMV - prEgivenMR)) / prEgivenMR;
		return Math.log10(likelihoodRatio);
	}
	private double llr(double prEgivenMR, double prEgivenMV, int mapq) {
		return llr(prEgivenMR, prEgivenMV, MathUtil.mapqToPr(mapq));
	}
	private double llr(double prEgivenMR, double prEgivenMV, int mapq1, int mapq2) {
		return llr(prEgivenMR, prEgivenMV, MathUtil.mapqToPr(mapq1, mapq2));
	}
	
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2) {
		double prEgivenMR = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength)); 
		double prEgivenMV = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, 0));
		return llr(prEgivenMR, prEgivenMV, mapq1, mapq2);
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq) {
		double prEgivenMR = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength)); 
		double prEgivenMV = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, 0));
		return llr(prEgivenMR, prEgivenMV, mapq);
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq) {
		double prEgivenMR = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(op, length)); 
		double prEgivenMV = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(op, 0));
		return llr(prEgivenMR, prEgivenMV, mapq);
	}
	
	public static double readPairFoldedCumulativeDistribution(IdsvSamFileMetrics metrics, int fragmentSize) {
		double pairsFromFragmentDistribution = 0;
		InsertSizeDistribution isd = metrics.getInsertSizeDistribution();
		IdsvMetrics im = metrics.getIdsvMetrics();
		if (fragmentSize > 0) {
			if (fragmentSize >= isd.getSupportLowerBound() && fragmentSize <= isd.getSupportUpperBound()) {
				double prUpper = 1.0 - isd.cumulativeProbability(fragmentSize - 1);
				double prLower = isd.cumulativeProbability(fragmentSize);
				double pr = Math.min(prUpper, prLower);
				pairsFromFragmentDistribution = pr * isd.getTotalMappedPairs();
			}
		}
		double totalPairs = im.READ_PAIRS_BOTH_MAPPED;
		double dpPairs = totalPairs - isd.getTotalMappedPairs() + pairsFromFragmentDistribution;
		return dpPairs / totalPairs;
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2) {
		double prEgivenMR = readPairFoldedCumulativeDistribution(metrics, fragmentSize);
		double prEgivenMV = 0.5; // TODO: actually calculate the inferred variant fragment size
		return llr(prEgivenMR, prEgivenMV, mapq1, mapq2);
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq) {
		IdsvMetrics im = metrics.getIdsvMetrics();
		// completely unmapped read pairs are excluded for consistency with sc and dp calculation
		double readPairs = im.READ_PAIRS - im.READ_PAIRS_ZERO_MAPPED;
		double prEgivenMR = im.READ_PAIRS_ONE_MAPPED / readPairs;
		// we assume that in our variant case, the read correctly maps across the breakpoint
		double prEgivenMV = im.READ_PAIRS_BOTH_MAPPED / readPairs;
		return llr(prEgivenMR, prEgivenMV, mapq);
	}
}
