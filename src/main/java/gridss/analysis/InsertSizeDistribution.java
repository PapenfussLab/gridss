package gridss.analysis;

import au.edu.wehi.idsv.util.CachedEnumeratedIntegerDistribution;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Log;
import picard.analysis.InsertSizeMetrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Set;

public class InsertSizeDistribution extends CachedEnumeratedIntegerDistribution {
	private static final Log log = Log.getInstance(InsertSizeDistribution.class);
	/**
	 * 
	 */
	private static final long serialVersionUID = -2332020573213102253L;
	private final long total;
	public long getTotalMappedPairs() {
		return total;
	}
	public static InsertSizeDistribution create(File insertSizeMetricsFile) {
		InsertSizeDistribution result = null;
		if (insertSizeMetricsFile != null && insertSizeMetricsFile.exists()) {
			FileReader reader = null;
			try {
				reader = new FileReader(insertSizeMetricsFile);
				MetricsFile<InsertSizeMetrics,Integer> file = new MetricsFile<InsertSizeMetrics, Integer>();
				file.read(reader);
				result = InsertSizeDistribution.create(file.getHistogram());
			} catch (FileNotFoundException e) {
				log.error("Missing insert size distribution for ", insertSizeMetricsFile);
			} finally {
				CloserUtil.close(reader);
			}
			if (result == null) {
				log.debug("Unable to extract insert size distribution from ", insertSizeMetricsFile);
			}
		}
		return result;
	}
	public static InsertSizeDistribution create(Histogram<Integer> insertSizeHistogram) {
		if (insertSizeHistogram == null) return null;
		int[] insertSize = new int[insertSizeHistogram.size()];
		double[] count = new double[insertSizeHistogram.size()];
		double total = insertSizeHistogram.getSumOfValues();
		int i = 0;
		Set<Integer> keys = insertSizeHistogram.keySet();
		for (Integer key : keys) {
			insertSize[i] = key;
			count[i] = insertSizeHistogram.get(key).getValue();
			i++;
		}
		return new InsertSizeDistribution(insertSize, count, (long)total);
	}
	public InsertSizeDistribution(int[] singletons, double[] readCounts) {
		this(singletons, readCounts, (long)Arrays.stream(readCounts).sum());
	}
	private InsertSizeDistribution(int[] singletons, double[] readCounts, long readTotal) {
		super(singletons, readCounts);
		this.total = readTotal;
	}
	public double descendingCumulativeProbability(int x) {
		return 1.0 - cumulativeProbability(x - 1);
	}
}
