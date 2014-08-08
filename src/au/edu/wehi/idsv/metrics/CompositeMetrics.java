package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.SamPairUtil.PairOrientation;

import java.util.Set;

import com.google.common.collect.Sets;

public class CompositeMetrics implements RelevantMetrics {
	private int maxFragmentSize = 0;
	private int maxReadLength = 0;
	private Set<PairOrientation> po = Sets.newHashSet();
	public CompositeMetrics(Iterable<RelevantMetrics> metrics) {
		for (RelevantMetrics m : metrics) {
			maxFragmentSize = Math.max(m.getMaxFragmentSize(), maxFragmentSize);
			maxReadLength =  Math.max(m.getMaxReadLength(), maxReadLength);
			po.add(m.getPairOrientation());
		}
	}
	@Override
	public double getMedianFragmentSize() {
		throw new IllegalStateException("Cannot determine mediam fragment size for composite metrics");
	}
	@Override
	public double getFragmentSizeStdDev() {
		throw new IllegalStateException("Cannot determine fragment std dev for composite metrics");
	}
	@Override
	public InsertSizeDistribution getInsertSizeDistribution() {
		throw new IllegalStateException("Cannot determine insert size distribution for composite metrics");
	}
	@Override
	public int getMaxFragmentSize() {
		return maxFragmentSize;
	}
	@Override
	public PairOrientation getPairOrientation() {
		if (po.size() != 1) throw new IllegalStateException("No consistent pair orientiation could be determined");
		return po.iterator().next();
	}
	@Override
	public int getMaxReadLength() {
		return maxReadLength;
	}
}
