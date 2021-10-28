package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.util.MathUtil;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import gridss.analysis.IdsvMetrics;
import htsjdk.samtools.CigarOperator;

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
	private double llr(double prEgivenMR, double prEgivenMV, int... mapq) {
		return llr(prEgivenMR, prEgivenMV, 1 - MathUtil.phredToPr(MathUtil.phredOr(Doubles.toArray(Ints.asList(mapq)))));
	}
	
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq1, int mapq2) {
		double prEgivenMR = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength)); 
		double prEgivenMV = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, 0));
		return llr(prEgivenMR, prEgivenMV, mapq1, mapq2);
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq) {
		double prEgivenMR = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, softclipLength)); 
		double prEgivenMV = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(CigarOperator.SOFT_CLIP, 0));
		return llr(prEgivenMR, prEgivenMV, mapq);
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, DirectedEvidence e, CigarOperator op, int length, int mapq) {
		double prEgivenMR = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(op, length)); 
		double prEgivenMV = MathUtil.phredToPr(metrics.getCigarDistribution().getPhred(op, 0));
		return llr(prEgivenMR, prEgivenMV, mapq);
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, DirectedEvidence e, int fragmentSize, int mapq1, int mapq2) {
		double prEgivenMR = MathUtil.phredToPr(metrics.getReadPairPhred(fragmentSize));
		double prEgivenMV = 0.5; // TODO: actually calculate the inferred variant fragment size
		return llr(prEgivenMR, prEgivenMV, mapq1, mapq2);
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, DirectedEvidence e, int mapq) {
		IdsvMetrics im = metrics.getIdsvMetrics();
		// completely unmapped read pairs are excluded for consistency with sc and dp calculation
		double readPairs = im.READ_PAIRS - im.READ_PAIRS_ZERO_MAPPED;
		double prEgivenMR = im.READ_PAIRS_ONE_MAPPED / readPairs;
		// we assume that in our variant case, the read correctly maps across the breakpoint
		double prEgivenMV = 0.5 * im.READ_PAIRS_BOTH_MAPPED / readPairs; // TODO: actually calculate the inferred variant fragment size
		return llr(prEgivenMR, prEgivenMV, mapq);
	}
}
