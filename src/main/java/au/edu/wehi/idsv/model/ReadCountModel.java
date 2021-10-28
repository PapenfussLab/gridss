package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
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
	public double scoreSplitRead(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq1, int mapq2) {
		return 1;
	}

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq) {
		return 1;
	}
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, DirectedEvidence e, CigarOperator op, int length, int mapq) {
		return 1;
	}

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, DirectedEvidence e, int fragmentSize, int mapq1, int mapq2) {
		return 1;
	}

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, DirectedEvidence e, int mapq) {
		return 1;
	}
	@Override
	public double scoreAssembly(DirectedEvidence e, int rp, double rpq, int sc, double scq, int localMapq, int remoteMapq) {
		return rp + sc;
	}
	@Override
	public double scoreBreakendAssembly(DirectedEvidence e, int rp, double rpq, int sc, double scq, int localMapq) {
		return rp + sc;
	}
}
