package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import htsjdk.samtools.CigarOperator;

public class ExclusionModel implements VariantScoringModel {
	private final VariantScoringModel model;
	private final boolean excludeDiscordantReadPairs;
	private final boolean excludeUnmappedMates;
	private final boolean excludeSplitReads;
	private final boolean excludeSoftClips;
	private final boolean excludeIndels;
	public ExclusionModel(
			VariantScoringModel model,
			boolean excludeDiscordantReadPairs,
			boolean excludeUnmappedMates,
			boolean excludeSplitReads,
			boolean excludeSoftClips,
			boolean excludeIndels) {
		this.model = model;
		this.excludeDiscordantReadPairs = excludeDiscordantReadPairs;
		this.excludeUnmappedMates = excludeUnmappedMates;
		this.excludeSplitReads = excludeSplitReads;
		this.excludeSoftClips = excludeSoftClips;
		this.excludeIndels = excludeIndels;
	}
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq1, int mapq2) {
		if (excludeSplitReads) return -1;
		return model.scoreSplitRead(metrics, e, softclipLength, mapq1, mapq2);
	}
	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq) {
		if (excludeSoftClips) return -1;
		return model.scoreSoftClip(metrics, e, softclipLength, mapq);
	}
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, DirectedEvidence e, CigarOperator op, int length, int mapq) {
		if (excludeIndels) return -1;
		return model.scoreIndel(metrics, e, op, length, mapq);
	}
	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, DirectedEvidence e, int fragmentSize, int mapq1, int mapq2) {
		if (excludeDiscordantReadPairs) return -1;
		return model.scoreReadPair(metrics, e, fragmentSize, mapq1, mapq2);
	}
	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics,DirectedEvidence e,  int mapq) {
		if (excludeUnmappedMates) return -1;
		return model.scoreUnmappedMate(metrics, e, mapq);
	}
}
