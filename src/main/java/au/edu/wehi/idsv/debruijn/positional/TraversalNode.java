package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;

import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public class TraversalNode {
	public final KmerPathSubnode node;
	public final int score;
	/**
	 * Length of path in kmers
	 */
	public final int pathLength;
	public final TraversalNode parent;
	/**
	 * Starting traversal node. Note that the root positional interval
	 * can be wider than that of the head.
	 */
	public final KmerPathSubnode root;
	public TraversalNode(KmerPathSubnode node, int baseScore) {
		this.node = node;
		this.score = baseScore + node.weight();
		this.parent = null;
		this.pathLength = node.length();
		this.root = node;
	}
	public TraversalNode(TraversalNode prev, KmerPathSubnode node) {
		this.node = node;
		this.score = prev.score + node.weight();
		this.parent = prev;
		this.pathLength = prev.pathLength + node.length();
		this.root = prev.root;
	}
	public TraversalNode(TraversalNode prev, KmerPathSubnode node, int terminalScore) {
		this.node = node;
		this.score = prev.score + node.weight() + terminalScore;
		this.parent = prev;
		this.pathLength = prev.pathLength + node.length();
		this.root = prev.root;
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
		this.root = node.root;
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
	public static Ordering<TraversalNode> ByFirstStart = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return Ints.compare(left.node.firstStart(), right.node.firstStart());
		}
	};
	public static Ordering<TraversalNode> ByLastEnd = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return Ints.compare(left.node.lastEnd(), right.node.lastEnd());
		}
	};
	public static Ordering<TraversalNode> ByScoreDescPosition = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(right.score, left.score)
					.compare(left.node.lastStart(), right.node.lastStart())
					.result();
		}
	};
}