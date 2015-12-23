package au.edu.wehi.idsv.model;

import htsjdk.samtools.CigarOperator;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;

/**
 * High-speed variant scoring model. This model assumes that one of the terms will dominate
 * thus the remainder can be ignored.
 * @author Daniel Cameron
 *
 */
public class MapqModel implements VariantScoringModel {
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2) {
		return MathUtil.prToPhred(1 - MathUtil.mapqToPr(mapq1, mapq2));
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq) {
		return MathUtil.prToPhred(1 - MathUtil.mapqToPr(mapq));
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq) {
		return MathUtil.prToPhred(1 - MathUtil.mapqToPr(mapq));
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2) {
		return MathUtil.prToPhred(1 - MathUtil.mapqToPr(mapq1, mapq2));
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq) {
		return MathUtil.prToPhred(1 - MathUtil.mapqToPr(mapq));
	}
}
