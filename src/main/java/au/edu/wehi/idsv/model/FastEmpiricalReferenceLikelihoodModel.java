package au.edu.wehi.idsv.model;

import htsjdk.samtools.CigarOperator;
import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.util.MathUtil;

/**
 * High-speed variant scoring model. This model assumes that one of the terms will dominate
 * thus the remainder can be ignored.
 * @author Daniel Cameron
 *
 */
public class FastEmpiricalReferenceLikelihoodModel implements VariantScoringModel {
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2) {
		double score = metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength);
		score = Math.min(score, mapq1);
		score = Math.min(score, mapq2);
		return score;
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq) {
		double score = metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength);
		score = Math.min(score, mapq);
		return score;
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq) {
		double score = metrics.getCigarDistribution().getPhred(op, length);
		score = Math.min(score, mapq);
		// score for indel is split over both sides of the event 
		return score;
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2) {
		double score = MathUtil.prToPhred(metrics.getReadPairPhred(fragmentSize));
		score = Math.min(score, mapq1);
		score = Math.min(score, mapq2);
		return score;
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq) {
		IdsvMetrics im = metrics.getIdsvMetrics();
		// completely unmapped read pairs are excluded for consistency with sc and dp calculation
		double score = MathUtil.prToPhred((double)im.READ_PAIRS_ONE_MAPPED / (double)(im.READ_PAIRS - im.READ_PAIRS_ZERO_MAPPED));
		score = Math.min(score, mapq);
		return score;
	}
}
