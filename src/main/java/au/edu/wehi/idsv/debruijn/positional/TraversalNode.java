package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class TraversalNode {
	public final KmerPathSubnode node;
	public final int score;
	public final int pathLength;
	public final TraversalNode parent;
	public TraversalNode(KmerPathSubnode node, int baseScore) {
		this.node = node;
		this.score = baseScore + node.weight();
		this.parent = null;
		this.pathLength = node.length();
	}
	public TraversalNode(TraversalNode prev, KmerPathSubnode node) {
		this.node = node;
		this.score = prev.score + node.weight();
		this.parent = prev;
		this.pathLength = prev.pathLength + node.length();
	}
	@Override
	public String toString() {
		String s = String.format("Score %d, %s", score, node);
		int length = getNodeCount();
		if (length > 2) {
			s = s + String.format(" -> (%d)", length);
		}
		if (length > 1) {
			s = s + getRoot().node().toString();
		}
		return s;
	}
	private KmerPathSubnode getRoot() {
		TraversalNode n = this;
		while (n.parent != null) n = n.parent;
		return n.node;
	}
	private int getNodeCount() {
		int count = 1;
		TraversalNode n = this;
		while (n.parent != null) {
			n = n.parent;
			count++;
		}
		return count;
	}
	/**
	 * Restricts the interval of an existing node to the given interval
	 * @param node existing node
	 * @param start new first kmer start position
	 * @param end new first kmer end position
	 */
	public TraversalNode(TraversalNode node, int start, int end) {
		assert(end >= start);
		assert(node.node.firstStart() <= start);
		assert(node.node.firstEnd() >= end);
		this.node = new KmerPathSubnode(node.node.node(), start, end);
		this.parent = node.parent;
		this.score = node.score;
		this.pathLength = node.pathLength;
	}
	/**
	 * Converts to a Subnode path under the assumption that the traversal was
	 * to successor nodes
	 * @return
	 */
	public ArrayDeque<KmerPathSubnode> toSubnodeNextPath() {
		ArrayDeque<KmerPathSubnode> contigPath = new ArrayDeque<KmerPathSubnode>();
		KmerPathSubnode last = node;
		contigPath.addFirst(node);
		for (TraversalNode n = parent; n != null; n = n.parent) {
			KmerPathSubnode current = n.node.givenNext(last);
			last = current;
			contigPath.addFirst(current);
		}
		return contigPath;
	}
	/**
	 * Converts to a Subnode path under the assumption that the traversal was
	 * to predecessor nodes
	 * @return
	 */
	public ArrayDeque<KmerPathSubnode> toSubnodePrevPath() {
		ArrayDeque<KmerPathSubnode> contigPath = new ArrayDeque<KmerPathSubnode>();
		KmerPathSubnode last = node;
		contigPath.addLast(node);
		for (TraversalNode n = parent; n != null; n = n.parent) {
			KmerPathSubnode current = n.node.givenPrev(last);
			last = current;
			contigPath.addLast(current);
		}
		return contigPath;
	}
	public static Ordering<TraversalNode> ByFirstStart = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return Ints.compare(left.node.firstStart(), right.node.firstStart());
		}
	};
	public  static Ordering<TraversalNode> ByLastEnd = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return Ints.compare(left.node.lastEnd(), right.node.lastEnd());
		}
	};
}