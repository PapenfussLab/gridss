package au.edu.wehi.socrates;

import java.io.File;

import net.sf.picard.analysis.CollectInsertSizeMetrics;
import net.sf.picard.analysis.InsertSizeMetrics;
import net.sf.picard.analysis.MetricAccumulationLevel;
import net.sf.picard.analysis.directed.InsertSizeMetricsCollector;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CollectionUtil;

public class RelevantMetrics {
	private InsertSizeMetrics insertSize = null;
	/**
	 * Creates a metric collector to record metrics required by Socrates
	 * @param header SAM header of file to process
	 * @return metric collector
	 */
	public static InsertSizeMetricsCollector createCollector(
			SAMFileHeader header) {
		return new InsertSizeMetricsCollector(
    			CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS, MetricAccumulationLevel.SAMPLE),
				header.getReadGroups(),
				// match CollectInsertSizeMetrics defaults
				 new CollectInsertSizeMetrics().MINIMUM_PCT,
				 new CollectInsertSizeMetrics().HISTOGRAM_WIDTH,
				 new CollectInsertSizeMetrics().DEVIATIONS);
	}
	public static void save(
			InsertSizeMetricsCollector metrics,
			MetricsFile<InsertSizeMetrics, Integer> metricsFile,
			File file) {		
		metrics.addAllLevelsToFile(metricsFile);
		metricsFile.write(file);
	}
	public RelevantMetrics(File file) {
		for (MetricBase metric : MetricsFile.readBeans(file)) {
			if (metric.getClass() == InsertSizeMetrics.class) {
				InsertSizeMetrics m = (InsertSizeMetrics)metric;
				if (m.SAMPLE == null &&
					m.LIBRARY == null &&
					m.READ_GROUP == null) {
					insertSize = m;
				}
			}
		}
	}
	public double getMedianFragmentSize() {
		// TODO: is this 5' difference or frag size?
		return insertSize.MEDIAN_INSERT_SIZE;
	}
	public double getFragmentSizeStdDev() {
		return 1.4826 * insertSize.MEDIAN_ABSOLUTE_DEVIATION;
	}
}
