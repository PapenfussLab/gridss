package au.edu.wehi.idsv.debruijn.positional;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

/**
 * KmerNode representing an individual kmer of a KmerPathNode
 * @author cameron.d
 *
 */
public class KmerPathNodeKmerNode implements KmerNode {
	private final KmerPathNode node;
	/**
	 * Offset indicates which kmer in the KmerPathNode this kmer corresponds to
	 * 
	 * A negative value indicates a collapsed kmer
	 */
	private final int offset;
	public KmerPathNodeKmerNode(KmerPathNode node, int offset) {
		this.node = node;
		this.offset = offset;
		assert(offset < node.length());
	}
	public KmerPathNodeKmerNode(int alternateKmerOffset, KmerPathNode node) {
		this.node = node;
		this.offset = -alternateKmerOffset - 1;
		assert(offsetOfPrimaryKmer() < node.length());
	}
	private int alternateKmerIndex() {
		return -offset - 1;
	}
	public KmerPathNode node() {
		return node;
	}
	public int offsetOfPrimaryKmer() {
		if (offset >= 0) return offset;
		return node.collapsedKmerOffsets().getInt(alternateKmerIndex());
	}
	@Override
	public long kmer() {
		if (offset >= 0) return node.kmer(offset);
		return node.collapsedKmers().getLong(alternateKmerIndex());
	}
	@Override
	public int startPosition() {
		return node.startPosition(offsetOfPrimaryKmer());
	}
	@Override
	public int endPosition() {
		return node.endPosition(offsetOfPrimaryKmer());
	}
	@Override
	public int weight() {
		return node.weight(offsetOfPrimaryKmer());
	}
	@Override
	public boolean isReference() {
		return node.isReference();
	}
	@Override
	public String toString() {
		return String.format("{%d}%s", offset, node);
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + node.hashCode();
		result = prime * result + offset;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		KmerPathNodeKmerNode other = (KmerPathNodeKmerNode) obj;
		if (other == null)
			return false;
		return offset == other.offset && node == other.node; 
	}
	/**
	 * Order by KmerPathNode, then by offset within the KmerPathNode
	 */
	public static Ordering<KmerPathNodeKmerNode> ByKmerPathNodeOffset = new Ordering<KmerPathNodeKmerNode>() {
		@Override
		public int compare(KmerPathNodeKmerNode left, KmerPathNodeKmerNode right) {
			return ComparisonChain.start()
					.compare(left.node, right.node, KmerNodeUtil.ByEndStartKmerReference)
					.compare(left.offsetOfPrimaryKmer(), right.offsetOfPrimaryKmer())
					.result();
		}
	};
}
