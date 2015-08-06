package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import picard.analysis.InsertSizeMetrics;
import au.edu.wehi.idsv.ProcessingContext;

import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

public class IdsvSamFileMetrics {
	private static final Log log = Log.getInstance(IdsvSamFileMetrics.class);
	public IdsvSamFileMetrics(ProcessingContext pc, File source) {
		this(pc.getFileSystemContext().getInsertSizeMetrics(source), pc.getFileSystemContext().getIdsvMetrics(source), pc.getFileSystemContext().getSoftClipMetrics(source));
	}
	public IdsvSamFileMetrics(File insertSizeMetricsFile, File idsvMetricsFile, File softClipMetricsFile) {
		this(getInsertSizeMetrics(insertSizeMetricsFile), getIdsvMetrics(idsvMetricsFile), getInsertSizeDistribution(insertSizeMetricsFile), getSoftClipMetrics(softClipMetricsFile));
	}
	private InsertSizeMetrics insertSize = null;
	private IdsvMetrics idsvMetrics = null;
	private InsertSizeDistribution insertDistribution = null;
	private List<SoftClipDetailMetrics> softClipDetailMetrics = null;
	private SoftClipSizeDistribution softClipDistribution;
	public IdsvMetrics getIdsvMetrics() { return idsvMetrics; }
	public InsertSizeMetrics getInsertSizeMetrics() { return insertSize; }
	public List<SoftClipDetailMetrics> getSoftClipDetailMetrics() { return softClipDetailMetrics; }
	
	private static InsertSizeDistribution getInsertSizeDistribution(File insertSizeMetricsFile) {
		return InsertSizeDistribution.create(insertSizeMetricsFile);
	}
	private static IdsvMetrics getIdsvMetrics(File idsvMetricsFiles) {
		IdsvMetrics metric = Iterators.getOnlyElement(Iterables.filter(MetricsFile.readBeans(idsvMetricsFiles), IdsvMetrics.class).iterator(), null);
		if (metric == null) {
			metric = new IdsvMetrics();
			log.error(String.format("No idsv metrics found in %s.", idsvMetricsFiles));
		}
		return metric;
	}
	public static InsertSizeMetrics getInsertSizeMetrics(File insertSizeMetricsFile) {
		InsertSizeMetrics insertSize = null;
		for (InsertSizeMetrics metric : Iterables.filter(MetricsFile.readBeans(insertSizeMetricsFile), InsertSizeMetrics.class)) {
			if (metric.SAMPLE == null && metric.LIBRARY == null && metric.READ_GROUP == null) {
				if (insertSize == null || insertSize.READ_PAIRS < metric.READ_PAIRS) {
					insertSize = metric;
				}
			}
		}
		if (insertSize == null) {
			insertSize = new InsertSizeMetrics();
			log.error(String.format("No pair-end insert size metrics found in %s.", insertSizeMetricsFile));
		}
		return insertSize;
	}
	public IdsvSamFileMetrics(InsertSizeMetrics insertSize, IdsvMetrics idsvMetrics, InsertSizeDistribution insertDistribution, List<SoftClipDetailMetrics> softClipDetailMetrics) {
		this.insertSize = insertSize;
		this.idsvMetrics = idsvMetrics;
		this.insertDistribution = insertDistribution;
		this.softClipDetailMetrics = softClipDetailMetrics;
		this.softClipDistribution = new SoftClipSizeDistribution(softClipDetailMetrics);
	}
	private static List<SoftClipDetailMetrics> getSoftClipMetrics(File softClipMetricsFile) {
		List<SoftClipDetailMetrics> sc = new ArrayList<SoftClipDetailMetrics>();
		for (SoftClipDetailMetrics metric : Iterables.filter(MetricsFile.readBeans(softClipMetricsFile), SoftClipDetailMetrics.class)) {
			sc.add(metric);
		}
		Collections.sort(sc, SoftClipDetailMetricsByLength);
		for (int i = 0; i < sc.size(); i++) {
			if (sc.get(i).LENGTH != i) {
				String msg = String.format("%s missing soft clip metric: expected metric for length %d, found %d", softClipMetricsFile, i, sc.get(0).LENGTH);
				log.error(msg);
				throw new RuntimeException(msg);
			}
		}
		return sc;
	}
	public InsertSizeDistribution getInsertSizeDistribution() {
		return insertDistribution;
	}
	public SoftClipSizeDistribution getSoftClipDistribution() {
		return softClipDistribution;
	}
	/**
	 * Sort order of soft clips by soft clip length
	 */
	private static Ordering<SoftClipDetailMetrics> SoftClipDetailMetricsByLength = new Ordering<SoftClipDetailMetrics>() {
		@Override
		public int compare(SoftClipDetailMetrics left, SoftClipDetailMetrics right) {
			return Longs.compare(left.LENGTH, right.LENGTH);
		}
	};
}
