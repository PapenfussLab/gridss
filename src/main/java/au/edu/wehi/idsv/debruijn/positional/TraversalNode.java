package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;

import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

public class TraversalNode {
	public final KmerPathSubnode node;
	public final int score;
	/**
	 * Length of path in kmers
	 */
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
	public TraversalNode(TraversalNode prev, KmerPathSubnode node, int additionalScore) {
		this.node = node;
		this.score = prev.score + node.weight() + additionalScore;
		this.parent = prev;
		this.pathLength = prev.pathLength + node.length();
	}
	/**
	 * First starting position of the first kmer in the path
	 * @return first position
	 */
	public int pathFirstStart() {
		return node.firstStart() - pathLength; 
	}
	/**
	 * Number of nodes in this path
	 * @return
	 */
	public int nodeCount() {
		int count = 0;
		TraversalNode n = this;
		while (n != null) {
			n = n.parent;
			count++;
		}
		return count;
	}
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(String.format("score=%d pathlength=%d\n", score, pathLength));
		for (TraversalNode n = this; n != null; n = n.parent) {
			sb.append(n.node.toString());
			//sb.append('\n');
		}
		return sb.toString();
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
	 * Determines whether this path already traverses through the given node
	 * @param node node to check traversal
	 * @return true if the path traverses this found, false otherwise
	 */
	public boolean traversingWouldCauseSelfIntersection(KmerPathNode node) {
		for (TraversalNode n = this; n != null; n = n.parent) {
			KmerPathNode pn = n.node.node();
			if (pn == node) return true;
			if (!IntervalUtil.overlapsClosed(pn.firstStart(), pn.firstEnd(), node.firstStart(), node.firstEnd())) {
				return false;
			}
		}
		return false;
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
	public static Ordering<TraversalNode> ByFirstStart = KmerNodeUtil.ByFirstStart.onResultOf((TraversalNode tn) -> tn.node);
	public static Ordering<TraversalNode> ByLastEndKmer = KmerNodeUtil.ByLastEndKmer.onResultOf((TraversalNode tn) -> tn.node);
	public static Ordering<TraversalNode> ByPathFirstStartArbitrary = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(left.pathFirstStart(), right.pathFirstStart())
					.compare(left, right, Ordering.arbitrary())
					.result();
		}
	};
	public static Ordering<TraversalNode> ByScoreDescPathFirstStartArbitrary = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(right.score, left.score)
					.compare(left, right, ByPathFirstStartArbitrary)
					.result();
		}
	};
}