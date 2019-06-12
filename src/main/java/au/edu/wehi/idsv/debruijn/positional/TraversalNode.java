package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.SanityCheckFailureException;
import au.edu.wehi.idsv.util.IntervalUtil;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

import java.util.ArrayDeque;

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
		return node.firstStart() - pathLength + node.length();
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
		sb.append(this.node.toString());
		if (this.parent != null) {
			sb.append(this.parent.node.toString());
			if (this.parent.parent != null) {
				TraversalNode n = this.parent.parent;
				while (n.parent != null) n = n.parent;
				sb.append("...\n");
				sb.append(n.node.toString());
			}
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
		checkValid(this);
		ArrayDeque<KmerPathSubnode> contigPath = new ArrayDeque<KmerPathSubnode>();
		KmerPathSubnode last = node;
		contigPath.addFirst(node);
		for (TraversalNode n = parent; n != null; n = n.parent) {
			checkValid(n);
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
		checkValid(this);
		ArrayDeque<KmerPathSubnode> contigPath = new ArrayDeque<KmerPathSubnode>();
		KmerPathSubnode last = node;
		contigPath.addLast(node);
		for (TraversalNode n = parent; n != null; n = n.parent) {
			checkValid(n);
			KmerPathSubnode current = n.node.givenPrev(last);
			last = current;
			contigPath.addLast(current);
		}
		return contigPath;
	}
	private void checkValid(TraversalNode tn) {
		if (!tn.node.node().isValid()) {
			throw new SanityCheckFailureException(String.format("Traversal of kmer length %d ending at [%d,%d] contains an invalid subnode over [%d-%d].",
					pathLength,
					node.lastStart(), node.lastEnd(),
					tn.node.firstStart(), tn.node.firstEnd()
					));
		}
	}
	public boolean sanityCheck() {
		assert(node != null);
		assert(node.sanityCheck());
		assert(score > 0);
		assert(node.node().isValid());
		if (parent != null) {
			assert(parent.sanityCheck());
			assert(score > parent.score);
			assert(parent.node.lastStart() + 1 <= node.firstStart());
			assert(parent.node.lastEnd() + 1 >= node.firstEnd());
			assert(pathLength == parent.pathLength + node.length());
		} else {
			assert(pathLength == node.length());
		}
		return true;
	}
	public static Ordering<TraversalNode> ByFirstStart = KmerNodeUtil.ByFirstStart.onResultOf((TraversalNode tn) -> tn.node);
	public static Ordering<TraversalNode> ByLastEndKmer = KmerNodeUtil.ByLastEndKmer.onResultOf((TraversalNode tn) -> tn.node);
	public static Ordering<TraversalNode> ByPathFirstStartScoreEndSubnode = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(left.pathFirstStart(), right.pathFirstStart())
					// JVisualVM indicates ArbitraryOrdering.identityHashCode() is a bottleneck
					.compare(left.score, right.score)
					.compare(left.node.firstStart(), right.node.firstStart())
					.compare(left.node.firstEnd(), right.node.firstEnd())
					.compare(left.node.firstKmer(), right.node.firstKmer())
					.result();
		}
	};
	public static Ordering<TraversalNode> ByScoreEndSubnode = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(left.score, right.score)
					.compare(left.node.firstStart(), right.node.firstStart())
					.compare(left.node.firstEnd(), right.node.firstEnd())
					.compare(left.node.firstKmer(), right.node.firstKmer())
					.result();
		}
	};
	public static Ordering<TraversalNode> ByKmerScoreStartEnd = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(left.node.firstKmer(), right.node.firstKmer())
					.compare(left.score, right.score)
					.compare(left.node.firstStart(), right.node.firstStart())
					.compare(left.node.firstEnd(), right.node.firstEnd())
					.result();
		}
	};
	public static Ordering<TraversalNode> ByScoreDescPathFirstEndSubnode = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(right.score, left.score)
					.compare(left, right, ByPathFirstStartScoreEndSubnode)
					.result();
		}
	};
}