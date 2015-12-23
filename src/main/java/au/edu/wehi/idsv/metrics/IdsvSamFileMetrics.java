package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import picard.analysis.InsertSizeMetrics;
import au.edu.wehi.idsv.ProcessingContext;

import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

public class IdsvSamFileMetrics {
	private static final Log log = Log.getInstance(IdsvSamFileMetrics.class);
	public IdsvSamFileMetrics(ProcessingContext pc, File source) {
		this(pc.getFileSystemContext().getInsertSizeMetrics(source), pc.getFileSystemContext().getIdsvMetrics(source), pc.getFileSystemContext().getCigarMetrics(source));
	}
	public IdsvSamFileMetrics(File insertSizeMetricsFile, File idsvMetricsFile, File cigarMetricsFile) {
		this(getInsertSizeMetrics(insertSizeMetricsFile), getIdsvMetrics(idsvMetricsFile), getInsertSizeDistribution(insertSizeMetricsFile), getCigarMetrics(cigarMetricsFile));
	}
	private InsertSizeMetrics insertSize = null;
	private IdsvMetrics idsvMetrics = null;
	private InsertSizeDistribution insertDistribution = null;
	private List<CigarDetailMetrics> cigarDetailMetrics = null;
	private CigarSizeDistribution cigarDistribution;
	public IdsvMetrics getIdsvMetrics() { return idsvMetrics; }
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
	public IdsvSamFileMetrics(InsertSizeMetrics insertSize, IdsvMetrics idsvMetrics, InsertSizeDistribution insertDistribution, List<CigarDetailMetrics> cigarDetailMetrics) {
		this.insertSize = insertSize;
		this.idsvMetrics = idsvMetrics;
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
	 * Sort order of soft clips by soft clip length
	 */
	private static Ordering<CigarDetailMetrics> CigarDetailMetricsByLength = new Ordering<CigarDetailMetrics>() {
		@Override
		public int compare(CigarDetailMetrics left, CigarDetailMetrics right) {
			return Longs.compare(left.LENGTH, right.LENGTH);
		}
	};
}
