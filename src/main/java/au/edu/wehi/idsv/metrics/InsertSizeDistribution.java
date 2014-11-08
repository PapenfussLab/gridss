package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Histogram.Bin;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Map.Entry;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;

import picard.analysis.InsertSizeMetrics;

public class InsertSizeDistribution extends EnumeratedIntegerDistribution {
	private static final Log log = Log.getInstance(InsertSizeDistribution.class);
	/**
	 * 
	 */
	private static final long serialVersionUID = -2332020573213102253L;
	private static final int MAX_CACHE_SIZE = 1024 * 1024;
	private double[] cpcache;
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
	public InsertSizeDistribution(int[] singletons, double[] probabilities, int readTotal) {
		super(singletons, probabilities);
		this.total = readTotal;
	}
	public double descendingCumulativeProbability(int x) {
		return 1.0 - cumulativeProbability(x - 1);
	}
	@Override
	public double cumulativeProbability(int x) {
		if (cpcache == null) {
			cpcache = new double[Math.min(MAX_CACHE_SIZE, getSupportUpperBound())];
			for (int i = 0; i < cpcache.length; i++) {
				cpcache[i] = super.cumulativeProbability(i);
			}
		}
		if (x < 0) return 0;
		if (x >= cpcache.length) return 1;
		return cpcache[x];
	}
}
