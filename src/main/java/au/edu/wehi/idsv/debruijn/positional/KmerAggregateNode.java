package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class KmerAggregateNode extends KmerNode {
	public long kmer() { return kmer; }
	public int startPosition() { return start; }
	public int endPosition() { return end; }
	public int weight() { return weight; }
	public boolean isReference() { return reference; }
	private final long kmer;
	private final int weight;
	private final int start;
	private final int end;
	private final boolean reference;
	public KmerAggregateNode(long kmer, int weight, int start, int end, boolean reference) {
		this.kmer = kmer;
		this.weight = weight;
		this.start = start;
		this.end = end;
		this.reference = reference;
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s%d, %s", isReference() ? "R" : " ",  startPosition(), endPosition(), weight(), KmerEncodingHelper.toApproximateString(kmer()));
	}
}
