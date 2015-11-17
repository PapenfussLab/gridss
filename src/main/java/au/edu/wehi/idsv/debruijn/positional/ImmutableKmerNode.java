package au.edu.wehi.idsv.debruijn.positional;

import java.util.Collection;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

/**
 * Total support for the given kmer over the given interval
 * @author Daniel Cameron
 *
 */
public class ImmutableKmerNode implements KmerNode {
	public long lastKmer() { return kmer; }
	public int lastStart() { return start; }
	public int lastEnd() { return end; }
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
		this(node.lastKmer(), node.lastEnd(), node.weight(), node.isReference(), node.lastStart());
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s%d, %s", lastStart(), lastEnd(), isReference() ? "R" : " ", weight(), KmerEncodingHelper.toApproximateString(lastKmer()));
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
		if (end != other.lastEnd())
			return false;
		if (kmer != other.lastKmer())
			return false;
		if (reference != other.isReference())
			return false;
		if (start != other.lastStart())
			return false;
		if (weight != other.weight())
			return false;
		return true;
	}
	/**
	 * Makes a copy of the given path node
	 * @param node
	 * @return
	 */
	public static Stream<ImmutableKmerNode> splitKmers(KmerNode node) {
		DeBruijnSequenceGraphNode sgn = (DeBruijnSequenceGraphNode)node;
		return IntStream.range(0, node.length())
				.mapToObj(i -> new ImmutableKmerNode(sgn.kmer(i), node.firstStart() + i, node.firstEnd() + i, node.isReference(), sgn.weight(i)));
	}
	/**
	 * Splits the given node into each individual position
	 * @param node
	 * @return
	 */
	public static Stream<ImmutableKmerNode> splitPositions(KmerNode node) {
		if (node instanceof DeBruijnSequenceGraphNode) throw new IllegalArgumentException("copyPath() should be called before fragment()");
		return IntStream.range(node.lastStart(), node.lastEnd() + 1)
				.mapToObj(p -> new ImmutableKmerNode(node.lastKmer(), p, p, node.isReference(), node.weight()));
	}
	/**
	 * Splits the given node into position de Bruijn graph nodes 
	 * @param list nodes
	 * @return individual constituent position de Bruijn graph nodes
	 */
	public static Stream<KmerNode> split(Collection<? extends KmerNode> list) {
		return list.stream().flatMap(n -> {
			if (n instanceof DeBruijnSequenceGraphNode) return ImmutableKmerNode.splitKmers(n).flatMap(x -> ImmutableKmerNode.splitPositions(x));
			return ImmutableKmerNode.splitPositions(n);
		});
	}
}
