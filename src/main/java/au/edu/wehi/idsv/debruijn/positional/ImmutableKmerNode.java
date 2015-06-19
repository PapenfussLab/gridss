package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class ImmutableKmerNode implements KmerNode {
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
	public ImmutableKmerNode(long kmer, int start, int end, boolean reference, int weight) {
		this.kmer = kmer;
		this.weight = weight;
		this.start = start;
		this.end = end;
		this.reference = reference;
	}
	public ImmutableKmerNode(KmerNode node) {
		this(node.kmer(), node.endPosition(), node.weight(), node.isReference(), node.startPosition());
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s%d, %s", startPosition(), endPosition(), isReference() ? "R" : " ", weight(), KmerEncodingHelper.toApproximateString(kmer()));
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + (int) (kmer ^ (kmer >>> 32));
		result = prime * result + (reference ? 1231 : 1237);
		result = prime * result + start;
		result = prime * result + weight;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		KmerNode other = (KmerNode) obj;
		if (end != other.endPosition())
			return false;
		if (kmer != other.kmer())
			return false;
		if (reference != other.isReference())
			return false;
		if (start != other.startPosition())
			return false;
		if (weight != other.weight())
			return false;
		return true;
	}
}
