package au.edu.wehi.idsv.model;

import htsjdk.samtools.CigarOperator;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;

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
	public double scoreSplitRead(IdsvSamFileMetrics metrics, int softclipLength, int mapq1, int mapq2) {
		if (excludeSplitReads) return -1;
		return model.scoreSplitRead(metrics, softclipLength, mapq1, mapq2);
	}
	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, int softclipLength, int mapq) {
		if (excludeSoftClips) return -1;
		return model.scoreSoftClip(metrics, softclipLength, mapq);
	}
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, CigarOperator op, int length, int mapq) {
		if (excludeIndels) return -1;
		return model.scoreIndel(metrics, op, length, mapq);
	}
	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, int fragmentSize, int mapq1, int mapq2) {
		if (excludeDiscordantReadPairs) return -1;
		return model.scoreReadPair(metrics, fragmentSize, mapq1, mapq2);
	}
	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, int mapq) {
		if (excludeUnmappedMates) return -1;
		return model.scoreUnmappedMate(metrics, mapq);
	}
}
