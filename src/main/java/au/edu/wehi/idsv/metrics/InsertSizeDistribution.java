package au.edu.wehi.idsv.metrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Map.Entry;

import au.edu.wehi.idsv.util.CachedEnumeratedIntegerDistribution;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Histogram.Bin;
import htsjdk.samtools.util.Log;
import picard.analysis.InsertSizeMetrics;

public class InsertSizeDistribution extends CachedEnumeratedIntegerDistribution {
	private static final Log log = Log.getInstance(InsertSizeDistribution.class);
	/**
	 * 
	 */
	private static final long serialVersionUID = -2332020573213102253L;
	private final int total;
	public int getTotalMappedPairs() {
		return total;
	}
	public static InsertSizeDistribution create(File insertSizeMetricsFile) {
		InsertSizeDistribution result = null;
		if (!insertSizeMetricsFile.exists()) return null;
		FileReader reader = null;
		try {
			reader = new FileReader(insertSizeMetricsFile);
			MetricsFile<InsertSizeMetrics,Integer> file = new MetricsFile<InsertSizeMetrics, Integer>();
			file.read(reader);
			result = InsertSizeDistribution.create(file.getHistogram());
		} catch (FileNotFoundException e) {
		} finally {
			CloserUtil.close(reader);
		}
		if (result == null) {
			log.error("Unable to extract insert size distribution from ", insertSizeMetricsFile);
		}
		return result;
	}
	public static InsertSizeDistribution create(Histogram<Integer> insertSizeHistogram) {
		if (insertSizeHistogram == null) return null;
		int[] insertSize = new int[insertSizeHistogram.size()];
		double[] count = new double[insertSizeHistogram.size()];
		double total = insertSizeHistogram.getSumOfValues();
		int i = 0;
		for (@SuppressWarnings("rawtypes") Entry<Integer, Bin> entry : insertSizeHistogram.entrySet()) {
			insertSize[i] = entry.getKey();
			count[i] = entry.getValue().getValue();// / total;
			i++;
		}
		return new InsertSizeDistribution(insertSize, count, (int)total);
	}
	public InsertSizeDistribution(int[] singletons, double[] readCounts, int readTotal) {
		super(singletons, readCounts);
		this.total = readTotal;
	}
	public double descendingCumulativeProbability(int x) {
		return 1.0 - cumulativeProbability(x - 1);
	}
}
