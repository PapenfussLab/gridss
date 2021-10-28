package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.util.MathUtil;
import gridss.analysis.IdsvMetrics;
import htsjdk.samtools.CigarOperator;

/**
 * Scoring model based on the likelihood of the evidence and
 * a correct mapping.
 * 
 * @author Daniel Cameron
 *
 */
public class EmpiricalReferenceLikelihoodModel implements VariantScoringModel {
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq1, int mapq2) {
		double score = MathUtil.phredOr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength), mapq1, mapq2);
		return score;
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq) {
		double score = MathUtil.phredOr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength), mapq);
		return score;
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, DirectedEvidence e, CigarOperator op, int length, int mapq) {
		double score = MathUtil.phredOr(metrics.getCigarDistribution().getPhred(op, length), mapq);
		return score;
	}
	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, DirectedEvidence e, int fragmentSize, int mapq1, int mapq2) {
		double score = MathUtil.phredOr(metrics.getReadPairPhred(fragmentSize), mapq1, mapq2);
		return score;
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, DirectedEvidence e, int mapq) {
		IdsvMetrics im = metrics.getIdsvMetrics();
		// completely unmapped read pairs are excluded for consistency with sc and dp calculation
		double prEgivenRM = (double)im.READ_PAIRS_ONE_MAPPED / (double)(im.READ_PAIRS - im.READ_PAIRS_ZERO_MAPPED);
		double score = MathUtil.phredOr(MathUtil.prToPhred(prEgivenRM), mapq);
		return score;
	}
}
