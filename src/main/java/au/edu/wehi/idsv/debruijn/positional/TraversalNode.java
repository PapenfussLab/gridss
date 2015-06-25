package au.edu.wehi.idsv.debruijn.positional;

public class TraversalNode {
	public TraversalNode(KmerPathSubnode node, int baseScore) {
		this.node = node;
		this.score = baseScore + node.weight();
		this.prev = null;
	}
	public TraversalNode(TraversalNode prev, KmerPathSubnode node) {
		this.node = node;
		this.score = prev.score + node.weight();
		this.prev = prev;
	}
	/**
	 * Restricts the interval of an existing node to the given interval
	 * @param node existing node
	 * @param start new first kmer start position
	 * @param end new first kmer end position
	 */
	public TraversalNode(TraversalNode node, int start, int end) {
		assert(node.node.firstKmerStartPosition() >= start);
		assert(node.node.firstKmerEndPosition() <= end);
		this.node = new KmerPathSubnode(node.node.node(), start, end);
		this.prev = node.prev;
		this.score = node.score;
	}
	public final KmerPathSubnode node;
	public final int score;
	public final TraversalNode prev;
}