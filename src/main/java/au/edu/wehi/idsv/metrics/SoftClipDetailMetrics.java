package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.metrics.MetricBase;

public class SoftClipDetailMetrics extends MetricBase {
	/**
	 * Length of soft clip
	 */
	public int LENGTH;
	/**
	 * Number of reads with a soft clip of this length
	 */
	public int READCOUNT;
	/**
	 * Number of reads with a soft clip that looks like adapter sequence
	 */
	//public int ADAPTERCOUNT;
}
