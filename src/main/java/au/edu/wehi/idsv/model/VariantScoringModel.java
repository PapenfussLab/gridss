package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import htsjdk.samtools.CigarOperator;

public interface VariantScoringModel {
	double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2);
	double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq);
	double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq);
	double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2);
	double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq);
}
