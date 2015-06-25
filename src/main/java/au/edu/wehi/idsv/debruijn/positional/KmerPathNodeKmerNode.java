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
	}
	public KmerPathNodeKmerNode(int alternateKmerOffset, KmerPathNode node) {
		this.node = node;
		this.offset = -alternateKmerOffset - 1;
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
	/**
	 * Order by KmerPathNode, then by offset within the KmerPathNode
	 */
	public static Ordering<KmerPathNodeKmerNode> ByKmerPathNodeOffset = new Ordering<KmerPathNodeKmerNode>() {
		@Override
		public int compare(KmerPathNodeKmerNode left, KmerPathNodeKmerNode right) {
			return ComparisonChain.start()
					.compare(left, right, KmerNode.ByEndStartKmerReference)
					.compare(left.offsetOfPrimaryKmer(), right.offsetOfPrimaryKmer())
					.result();
		}
	};
}
