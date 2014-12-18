package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.metrics.MetricBase;

public class SoftClipDetailMetrics extends MetricBase {
	/**
	 * Length of soft clip
	 */
	public long LENGTH;
	/**
	 * Number of reads with a soft clip of this length
	 */
	public long READCOUNT;
	/**
	 * Number of reads with a soft clip that looks like adapter sequence
	 */
	//public long ADAPTERCOUNT;
}
