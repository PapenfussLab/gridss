package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import htsjdk.samtools.CigarOperator;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Custom evidence weighting model. The evidence weighting is encoded in the "weight" named group of the read name regex.
 */
public class ReadNameWeightedModel implements VariantScoringModel {
	private static final String CAPTURE_GROUP_NAME = "weight";
	private static final double NO_REGEX_MATCH_SCORE = 1;
	private final Pattern pattern;

	public ReadNameWeightedModel(String regex) {
		this.pattern = Pattern.compile(regex);
	}
	private double readWeight(DirectedEvidence e, double defaultScore) {
		Matcher match = pattern.matcher(e.getUnderlyingSAMRecord().getReadName());
		if (!match.matches()) return defaultScore;
		String group = match.group(CAPTURE_GROUP_NAME);
		if (group == null) return defaultScore;
		return Double.parseDouble(group);
	}
	@Override
	public double scoreSplitRead(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq1, int mapq2) { return readWeight(e, NO_REGEX_MATCH_SCORE); }

	@Override
	public double scoreSoftClip(IdsvSamFileMetrics metrics, DirectedEvidence e, int softclipLength, int mapq){ return readWeight(e, NO_REGEX_MATCH_SCORE); }
	
	@Override
	public double scoreIndel(IdsvSamFileMetrics metrics, DirectedEvidence e, CigarOperator op, int length, int mapq) { return readWeight(e, NO_REGEX_MATCH_SCORE); }

	@Override
	public double scoreReadPair(IdsvSamFileMetrics metrics, DirectedEvidence e, int fragmentSize, int mapq1, int mapq2) { return readWeight(e, NO_REGEX_MATCH_SCORE); }

	@Override
	public double scoreUnmappedMate(IdsvSamFileMetrics metrics, DirectedEvidence e, int mapq) { return readWeight(e, NO_REGEX_MATCH_SCORE); }
	@Override
	public double scoreAssembly(DirectedEvidence e, int rp, double rpq, int sc, double scq, int localMapq, int remoteMapq) {
		return readWeight(e, rpq + scq);
	}
	@Override
	public double scoreBreakendAssembly(DirectedEvidence e, int rp, double rpq, int sc, double scq, int localMapq) {
		return readWeight(e, rpq + scq);
	}
}
