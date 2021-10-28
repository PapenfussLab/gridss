package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.util.MathUtil;
import htsjdk.samtools.CigarOperator;

/**
 * High-speed variant scoring model. This model assumes that one of the terms will dominate
 * thus the remainder can be ignored.
 * @author Daniel Cameron
 *
 */
public class MapqModel implements VariantScoringModel {
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq1, int mapq2) {
		return MathUtil.phredOr(mapq1, mapq2);
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq) {
		return MathUtil.phredOr(mapq);
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, DirectedEvidence e, CigarOperator op, int length, int mapq) {
		return MathUtil.phredOr(mapq);
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, DirectedEvidence e, int fragmentSize, int mapq1, int mapq2) {
		return MathUtil.phredOr(mapq1, mapq2);
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, DirectedEvidence e, int mapq) {
		return MathUtil.phredOr(mapq);
	}
}
