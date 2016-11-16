package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import htsjdk.samtools.CigarOperator;

/**
 * High-speed variant scoring model. This model assumes that one of the terms will dominate
 * thus the remainder can be ignored.
 * @author Daniel Cameron
 *
 */
public class ReadCountModel implements VariantScoringModel {
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2) {
		return 1;
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq) {
		return 1;
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq) {
		return 1;
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2) {
		return 1;
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq) {
		return 1;
	}
	@Override
	public double scoreAssembly(int rp, double rpq, int sc, double scq, int localMapq, int remoteMapq) {
		return rp + sc;
	}
	@Override
	public double scoreBreakendAssembly(int rp, double rpq, int sc, double scq, int localMapq) {
		return rp + sc;
	}
}
