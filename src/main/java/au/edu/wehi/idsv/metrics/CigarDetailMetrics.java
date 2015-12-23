package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.metrics.MetricBase;

public class CigarDetailMetrics extends MetricBase {
	/**
	 * Length of indel
	 */
	public int LENGTH;
	/**
	 * Number of reads with a indel of this length
	 */
	public long COUNT;
	/**
	 * Cigar operator
	 */
	public char OPERATOR;
}
