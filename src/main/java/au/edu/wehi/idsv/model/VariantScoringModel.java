package au.edu.wehi.idsv.model;

import htsjdk.samtools.CigarOperator;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;

public interface VariantScoringModel {
	double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2);
	double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq);
	double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq);
	double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2);
	double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq);
	default double scoreBreakendAssembly(int rp, double rpq, int sc, double scq, int localMapq) {
		double qual = rpq + scq;
		qual = Math.min(localMapq * (rp + sc), qual);
		return qual;
	}
	default double scoreAssembly(int rp, double rpq, int sc, double scq, int localMapq, int remoteMapq) {
		double qual = rpq + scq;
		qual = Math.min(localMapq * (rp + sc), qual);
		qual = Math.min(remoteMapq * (rp + sc), qual);
		return qual;
	}
}
