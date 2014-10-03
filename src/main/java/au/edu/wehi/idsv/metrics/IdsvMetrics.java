package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.metrics.MetricBase;

public class IdsvMetrics extends MetricBase {
	/**
	 * Length of longest read
	 */
	public Integer MAX_READ_LENGTH = 0;
	/**
	 * Inferred size of largest concordantly mapped fragment including all read bases of both pairs
	 * Note: This differs from the htsjdk definition of fragment size 
	 */
	public Integer MAX_PROPER_PAIR_FRAGMENT_LENGTH = 0; // FIXME: need to calc proper pairs ourselves as can't trust aligner
	public static IdsvMetrics merge(IdsvMetrics arg0, IdsvMetrics arg1) {
		IdsvMetrics m = new IdsvMetrics();
		m.MAX_READ_LENGTH = Math.max(arg0.MAX_READ_LENGTH, arg1.MAX_READ_LENGTH);
		m.MAX_PROPER_PAIR_FRAGMENT_LENGTH = Math.max(arg0.MAX_PROPER_PAIR_FRAGMENT_LENGTH, arg1.MAX_PROPER_PAIR_FRAGMENT_LENGTH);
		return m;
	}
}
