package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import htsjdk.samtools.CigarOperator;

public interface VariantScoringModel {
	double scoreSplitRead(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq1, int mapq2);
	double scoreSoftClip(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq);
	double scoreIndel(IdsvSamFileMetrics metrics, DirectedEvidence e, CigarOperator op, int length, int mapq);
	double scoreReadPair(IdsvSamFileMetrics metrics, DirectedEvidence e, int fragmentSize, int mapq1, int mapq2);
	double scoreUnmappedMate(IdsvSamFileMetrics metrics, DirectedEvidence e, int mapq);
	default double scoreBreakendAssembly(DirectedEvidence e, int rp, double rpq, int sc, double scq, int localMapq) {
		double qual = rpq + scq;
		qual = Math.min(localMapq * (rp + sc), qual);
		return qual;
	}
	default double scoreAssembly(DirectedEvidence e, int rp, double rpq, int sc, double scq, int localMapq, int remoteMapq) {
		double qual = rpq + scq;
		qual = Math.min(localMapq * (rp + sc), qual);
		qual = Math.min(remoteMapq * (rp + sc), qual);
		return qual;
	}
}
