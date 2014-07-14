package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.List;

import picard.analysis.CollectInsertSizeMetrics;
import picard.analysis.InsertSizeMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.directed.InsertSizeMetricsCollector;

import com.google.common.collect.ImmutableList;

public class PicardMetrics implements RelevantMetrics {
	private Log log = Log.getInstance(PicardMetrics.class);
	private InsertSizeMetrics insertSize = null;
	/**
	 * Creates a metric collector to record metrics required by idsv
	 * @param header SAM header of file to process
	 * @return metric collector
	 */
	public static InsertSizeMetricsCollector createCollector(
			SAMFileHeader header) {
		List<SAMReadGroupRecord> rg = ImmutableList.<SAMReadGroupRecord>of();
		if (header != null) {
			rg = header.getReadGroups();
		}
		return new InsertSizeMetricsCollector(
    			CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS), null, //, MetricAccumulationLevel.SAMPLE), rg,
				// match CollectInsertSizeMetrics defaults
				new CollectInsertSizeMetrics().MINIMUM_PCT,
				new CollectInsertSizeMetrics().Histogram_WIDTH,
				new CollectInsertSizeMetrics().DEVIATIONS);
	}
	public static void save(
			InsertSizeMetricsCollector metrics,
			MetricsFile<InsertSizeMetrics, Integer> metricsFile,
			File file) {
		metrics.finish();
		metrics.addAllLevelsToFile(metricsFile);
		metricsFile.write(file);
	}
	public PicardMetrics(File file) {
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
		if (insertSize == null) {
			insertSize = new InsertSizeMetrics();
			log.error(String.format("No pair-end insert size metrics found in %s.", file));
		}
	}
	protected PicardMetrics() {
		insertSize = new InsertSizeMetrics();
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.RelevantMetrics#getMedianFragmentSize()
	 */
	@Override
	public double getMedianFragmentSize() {
		// TODO: is this 5' difference or frag size?
		return insertSize.MEDIAN_INSERT_SIZE;
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.RelevantMetrics#getFragmentSizeStdDev()
	 */
	@Override
	public double getFragmentSizeStdDev() {
		return 1.4826 * insertSize.MEDIAN_ABSOLUTE_DEVIATION;
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.RelevantMetrics#getMaxFragmentSize()
	 */
	@Override
	public int getMaxFragmentSize() {
		// TODO: is this 5' difference or frag size?
		return (int)Math.ceil(getMedianFragmentSize() + 3 * getFragmentSizeStdDev());
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.RelevantMetrics#getPairOrientation()
	 */
	@Override
	public PairOrientation getPairOrientation() {
		if (insertSize == null) return null;
		return insertSize.PAIR_ORIENTATION;
	}
}
