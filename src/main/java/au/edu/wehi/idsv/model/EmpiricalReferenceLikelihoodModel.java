package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import htsjdk.samtools.CigarOperator;

/**
 * High-speed variant scoring model. This model assumes that one of the terms will dominate
 * thus the remainder can be ignored.
 * @author Daniel Cameron
 *
 */
public class EmpiricalReferenceLikelihoodModel implements VariantScoringModel {
	public double referenceLikelihood(double prEgivenRM, double prM) {
		double prEgivenRMbar = 1;
		return prM * prEgivenRM; // + (1 - prM) * prEgivenRMbar;
	}
	public double score(double prEgivenRM, int mapq) {
		return MathUtil.prToPhred(1 - referenceLikelihood(prEgivenRM, MathUtil.mapqToPr(mapq)));
	}
	public double score(double prEgivenRM, int mapq1, int mapq2) {
		return MathUtil.prToPhred(1 - referenceLikelihood(prEgivenRM, MathUtil.mapqToPr(mapq1, mapq2)));
	}
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2) {
		return score(MathUtil.mapqToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength)), mapq1, mapq2);
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq) {
		return score(MathUtil.mapqToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength)), mapq);
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq) {
		return score(MathUtil.mapqToPr(metrics.getCigarDistribution().getPhred(op, length)), mapq);
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2) {
		return score(EmpiricalLlrModel.readPairFoldedCumulativeDistribution(metrics, fragmentSize), mapq1, mapq2);
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq) {
		IdsvMetrics im = metrics.getIdsvMetrics();
		// completely unmapped read pairs are excluded for consistency with sc and dp calculation
		double prEgivenRM = (double)im.READ_PAIRS_ONE_MAPPED / (double)(im.READ_PAIRS - im.READ_PAIRS_ZERO_MAPPED);
		return score(prEgivenRM, mapq);
	}
}
