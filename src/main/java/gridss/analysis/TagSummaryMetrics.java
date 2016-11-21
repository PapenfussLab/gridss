package gridss.analysis;

import htsjdk.samtools.metrics.MetricBase;

public class TagSummaryMetrics extends MetricBase {
	/**
	 * tag name 
	 */
	public String TAG;
	/**
	 * Number of reads containing attribute
	 */
	public long COUNT;
}
