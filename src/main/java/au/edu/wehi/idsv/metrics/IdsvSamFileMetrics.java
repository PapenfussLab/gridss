package au.edu.wehi.idsv.metrics;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.util.MathUtil;
import gridss.analysis.CigarDetailMetrics;
import gridss.analysis.CigarSizeDistribution;
import gridss.analysis.IdsvMetrics;
import gridss.analysis.InsertSizeDistribution;
import gridss.analysis.MapqMetrics;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Log;
import picard.analysis.InsertSizeMetrics;

public class IdsvSamFileMetrics {
	private static final Log log = Log.getInstance(IdsvSamFileMetrics.class);
	public IdsvSamFileMetrics(GenomicProcessingContext pc, File source) {
		this(pc.getFileSystemContext().getInsertSizeMetrics(source),
				pc.getFileSystemContext().getIdsvMetrics(source),
				pc.getFileSystemContext().getMapqMetrics(source),
				pc.getFileSystemContext().getCigarMetrics(source));
	}
	public IdsvSamFileMetrics(File insertSizeMetricsFile, File idsvMetricsFile, File mapqMetricsFile, File cigarMetricsFile) {
		this(getInsertSizeMetrics(insertSizeMetricsFile),
				getIdsvMetrics(idsvMetricsFile),
				getMapqMetrics(mapqMetricsFile),
				getInsertSizeDistribution(insertSizeMetricsFile),
				getCigarMetrics(cigarMetricsFile));
	}
	private InsertSizeMetrics insertSize = null;
	private MapqMetrics mapqMetrics = null;
	private IdsvMetrics idsvMetrics = null;
	private InsertSizeDistribution insertDistribution = null;
	private List<CigarDetailMetrics> cigarDetailMetrics = null;
	private CigarSizeDistribution cigarDistribution;
	public IdsvMetrics getIdsvMetrics() { return idsvMetrics; }
	public MapqMetrics getMapqMetrics() { return mapqMetrics; }
	public InsertSizeMetrics getInsertSizeMetrics() { return insertSize; }
	public List<CigarDetailMetrics> getCigarDetailMetrics() { return cigarDetailMetrics; }
	
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
		InsertSizeMetrics bestMetrics = null;
		if (insertSizeMetricsFile.exists()) {
			for (InsertSizeMetrics metric : Iterables.filter(MetricsFile.readBeans(insertSizeMetricsFile), InsertSizeMetrics.class)) {
				if (metric.SAMPLE == null && metric.LIBRARY == null && metric.READ_GROUP == null) {
					if (bestMetrics == null || bestMetrics.READ_PAIRS < metric.READ_PAIRS) {
						bestMetrics = metric;
					}
				}
			}
		}
		if (bestMetrics == null) {
			log.warn(String.format("No pair-end insert size metrics found in %s. Assuming this library contains single-end reads", insertSizeMetricsFile));
		}
		return bestMetrics;
	}
	public static MapqMetrics getMapqMetrics(File mapqMetricsFile) {
		MapqMetrics bestMetrics = null;
		for (MapqMetrics metric : Iterables.filter(MetricsFile.readBeans(mapqMetricsFile), MapqMetrics.class)) {
			if (metric.SAMPLE == null && metric.LIBRARY == null && metric.READ_GROUP == null) {
				if (bestMetrics == null || bestMetrics.MAPPED_READS < metric.MAPPED_READS) {
					bestMetrics = metric;
				}
			}
		}
		if (bestMetrics == null) {
			bestMetrics = new MapqMetrics();
			log.error(String.format("No mapq metrics found in %s.", mapqMetricsFile));
		}
		return bestMetrics;
	}
	public IdsvSamFileMetrics(InsertSizeMetrics insertSize, IdsvMetrics idsvMetrics, MapqMetrics mapqMetrics, InsertSizeDistribution insertDistribution, List<CigarDetailMetrics> cigarDetailMetrics) {
		this.insertSize = insertSize;
		this.idsvMetrics = idsvMetrics;
		this.mapqMetrics = mapqMetrics;
		this.insertDistribution = insertDistribution;
		this.cigarDetailMetrics = cigarDetailMetrics;
		this.cigarDistribution = new CigarSizeDistribution(cigarDetailMetrics);
	}
	private static List<CigarDetailMetrics> getCigarMetrics(File cigarMetricsFile) {
		List<CigarDetailMetrics> list = new ArrayList<CigarDetailMetrics>();
		for (CigarDetailMetrics metric : Iterables.filter(MetricsFile.readBeans(cigarMetricsFile), CigarDetailMetrics.class)) {
			list.add(metric);
		}
		Collections.sort(list, CigarDetailMetricsByLength);
		for (CigarOperator op : CigarOperator.values()) {
			List<CigarDetailMetrics> opList = list.stream().filter(cdm -> cdm.OPERATOR == CigarOperator.enumToCharacter(op)).collect(Collectors.toList());
			for (int i = 0; i < opList.size(); i++) {
				if (opList.get(i).LENGTH != i) {
					String msg = String.format("%s missing cigar metric: expected metric for length %d, found %d", cigarMetricsFile, i, opList.get(0).LENGTH);
					log.error(msg);
					throw new RuntimeException(msg);
				}
			}
		}
		return list;
	}
	public InsertSizeDistribution getInsertSizeDistribution() {
		return insertDistribution;
	}
	public CigarSizeDistribution getCigarDistribution() {
		return cigarDistribution;
	}
	private int maxSoftClipLength = -1;
	public int getMaxSoftClipLength() {
		if (maxSoftClipLength < 0) {
			maxSoftClipLength = getCigarDetailMetrics().stream().filter(m -> m.OPERATOR == CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP)).mapToInt(m -> m.LENGTH).max().orElse(0);
		}
		return maxSoftClipLength;
	}
	/**
	 * Returns the phred-scaled likelihood of a fragment size at least as extreme as the given size.
	 * @param fragmentSize fragment size
	 * @return phred-scaled likelihood of a fragment as or more extreme
	 */
	public double getReadPairPhred(int fragmentSize) {
		return MathUtil.prToPhred(readPairFoldedCumulativeDistribution(fragmentSize));
	}
	public double readPairFoldedCumulativeDistribution(int fragmentSize) {
		double pairsFromFragmentDistribution = 0;
		if (fragmentSize > 0) {
			if (fragmentSize >= insertDistribution.getSupportLowerBound() && fragmentSize <= insertDistribution.getSupportUpperBound()) {
				double prUpper = 1.0 - insertDistribution.cumulativeProbability(fragmentSize - 1);
				double prLower = insertDistribution.cumulativeProbability(fragmentSize);
				double pr = Math.min(prUpper, prLower);
				pairsFromFragmentDistribution = pr * insertDistribution.getTotalMappedPairs();
			}
		}
		double totalPairs = idsvMetrics.READ_PAIRS_BOTH_MAPPED;
		double dpPairs = totalPairs - insertDistribution.getTotalMappedPairs() + pairsFromFragmentDistribution;
		return dpPairs / totalPairs;
	}
	
	/**
	 * Sort order of soft clips by soft clip length
	 */
	private static Ordering<CigarDetailMetrics> CigarDetailMetricsByLength = new Ordering<CigarDetailMetrics>() {
		@Override
		public int compare(CigarDetailMetrics left, CigarDetailMetrics right) {
			return Longs.compare(left.LENGTH, right.LENGTH);
		}
	};
}
