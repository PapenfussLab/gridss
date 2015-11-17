package au.edu.wehi.idsv.metrics;

import java.util.List;

/**
 * Distribution of soft clip lengths
 * @author Daniel Cameron
 *
 */
public class SoftClipSizeDistribution {
	private final double[] phred;
	public SoftClipSizeDistribution(List<SoftClipDetailMetrics> sc) {
		this.phred = calcPhred(sc);
	}
	private static double[] calcPhred(List<SoftClipDetailMetrics> sc) {
		long total = 0;
		for (SoftClipDetailMetrics m : sc) {
			total += m.READCOUNT;
		}
		double[] phred = new double[sc.size()];
		double score = 0;
		long cumsum = total;
		for (int i = 0; i < sc.size(); i++) {
			SoftClipDetailMetrics m = sc.get(i);
			assert(m.LENGTH == i);
			if (cumsum > 0) {
				score = -10 * Math.log10((double)cumsum / (double)total);
			}
			phred[i] = score;
			cumsum -= m.READCOUNT;
		}
		if (phred.length == 0) {
			phred = new double[] { 0 };
		}
		return phred;
	}
	/**
	 * Returns the phred scaled probability of a soft clip of at least this length
	 * @param softClipLength
	 * @return soft clip phred score
	 */
	public double getPhred(int softClipLength) {
		if (softClipLength < 0) return 0;
		if (softClipLength >= phred.length) return phred[phred.length - 1];
		return phred[softClipLength];
	}
}
