package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileReader;
import java.util.List;

import picard.analysis.CollectInsertSizeMetrics;
import picard.analysis.InsertSizeMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.directed.InsertSizeMetricsCollector;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;

public class IdsvSamFileMetrics implements RelevantMetrics {
	private Log log = Log.getInstance(IdsvSamFileMetrics.class);
	private InsertSizeMetrics insertSize = null;
	private IdsvMetrics idsvMetrics = null;
	private InsertSizeDistribution insertDistribution = null;
	/**
	 * Creates a metric collector to record metrics required by idsv
	 * @param header SAM header of file to process
	 * @return metric collector
	 */
	public static InsertSizeMetricsCollector createInsertSizeMetricsCollector(SAMFileHeader header) {
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
	/**
	 * Creates a metric collector to record metrics required by idsv
	 * @param header SAM header of file to process
	 * @return metric collector
	 */
	public static IdsvMetricsCollector createIdsvMetricsCollector() {
		return new IdsvMetricsCollector();
	}
	public IdsvSamFileMetrics(File insertSizeMetricsFile, File idsvMetricsFiles) {
		for (MetricBase metric : Iterables.concat(MetricsFile.readBeans(insertSizeMetricsFile), MetricsFile.readBeans(idsvMetricsFiles))) {
			if (metric.getClass() == InsertSizeMetrics.class) {
				InsertSizeMetrics m = (InsertSizeMetrics)metric;
				if (m.SAMPLE == null &&
					m.LIBRARY == null &&
					m.READ_GROUP == null) {
					insertSize = m;
				}
			} else if (metric.getClass() == IdsvMetrics.class) {
				idsvMetrics = (IdsvMetrics)metric;
			}
		}
		if (insertSize == null) {
			insertSize = new InsertSizeMetrics();
			log.error(String.format("No pair-end insert size metrics found in %s.", insertSizeMetricsFile));
		}
		if (idsvMetrics == null) {
			idsvMetrics = new IdsvMetrics();
			log.error(String.format("No idsv metrics found in %s.", idsvMetricsFiles));
		}
		insertDistribution = InsertSizeDistribution.create(insertSizeMetricsFile);
	}
	public IdsvSamFileMetrics(InsertSizeMetrics insertSizeMetrics, IdsvMetrics idsvMetrics) {
		if (insertSizeMetrics == null) throw new NullPointerException("insertSize");
		if (idsvMetrics == null) throw new NullPointerException("idsvMetrics");
		this.insertSize = insertSizeMetrics;
		this.idsvMetrics = idsvMetrics;
	}
	protected IdsvSamFileMetrics() {
		insertSize = new InsertSizeMetrics();
		idsvMetrics = new IdsvMetrics();
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
		// return (int)Math.ceil(getMedianFragmentSize() + 3 * getFragmentSizeStdDev());
		return Math.max(getMaxReadLength(), idsvMetrics.MAX_PROPER_PAIR_FRAGMENT_LENGTH);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.RelevantMetrics#getPairOrientation()
	 */
	@Override
	public PairOrientation getPairOrientation() {
		if (insertSize == null) return null;
		return insertSize.PAIR_ORIENTATION;
	}
	@Override
	public int getMaxReadLength() {
		return idsvMetrics.MAX_READ_LENGTH;
	}
	@Override
	public InsertSizeDistribution getInsertSizeDistribution() {
		return insertDistribution;
	}
}
