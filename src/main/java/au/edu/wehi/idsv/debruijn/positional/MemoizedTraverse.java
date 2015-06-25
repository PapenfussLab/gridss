package au.edu.wehi.idsv.debruijn.positional;

import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.TreeSet;

import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * Helper class to track memoization of nodes during positional graph traversal
 * @author Daniel Cameron
 *
 */
public class MemoizedTraverse {
	/**
	 * Since a positional de Bruijn graph is a directed acyclic graph,
	 * we can calculate maximal weighted paths by a positional traverse
	 * (BFS in position space) of the graph, caching the best predecessor
	 * of each node. 
	 */
	private final NavigableSet<TraversalNode> memoized = new TreeSet<TraversalNode>(MemoizedNodeByKmerStartEndScore);
	private final PriorityQueue<TraversalNode> frontier = new PriorityQueue<TraversalNode>(1024, MemoizedNodeByEnd);
	/**
	 * Memoize the given score for the given position
	 * @param score initial score
	 * @param start starting node
	 */
	public void memoize(TraversalNode node) {
		// check existing scores for the kmer interval
		NavigableSet<TraversalNode> couldOverlap = memoized;
		TraversalNode floorNode = memoized.floor(node);
		if (floorNode != null) {
			couldOverlap = memoized.tailSet(memoized.floor(node), true);
		}
		for (TraversalNode existing : couldOverlap) {
			if (existing.node.kmer(0) > node.node.kmer(0)) break;
			if (existing.node.kmer(0) < node.node.kmer(0)) continue;
			if (existing.node.firstKmerEndPosition() < node.node.firstKmerStartPosition()) continue;
			if (existing.node.firstKmerStartPosition() > node.node.firstKmerEndPosition()) break;
			
			// we overlap an existing path
			assert(existing.node.kmer(0) == node.node.kmer(0) && IntervalUtil.overlapsClosed(
					existing.node.firstKmerStartPosition(), existing.node.firstKmerStartPosition(),
					node.node.firstKmerStartPosition(), node.node.firstKmerStartPosition()));
			
			// ok, so now we know the nodes overlap
			if (node.score > existing.score) {
				// remove existing node in overlapping interval
				memoized.remove(existing);
				if (existing.node.firstKmerStartPosition() < node.node.firstKmerStartPosition()) {
					// still valid in earlier interval
					TraversalNode existingBefore = new TraversalNode(existing, existing.node.firstKmerStartPosition(), node.node.firstKmerStartPosition() - 1);
					memoized.add(existingBefore);
					frontier.add(existingBefore);
				}
				if (existing.node.firstKmerEndPosition() > node.node.firstKmerEndPosition()) {
					// still valid in earlier interval
					TraversalNode existingAfter = new TraversalNode(existing, node.node.firstKmerEndPosition() + 1, existing.node.firstKmerEndPosition());
					memoized.add(existingAfter);
					frontier.add(existingAfter);
				}
			} else { // existing node scores better than us
				int newStartPosition = existing.node.firstKmerEndPosition() + 1;
				if (node.node.firstKmerStartPosition() < existing.node.firstKmerStartPosition()) {
					// start before this node -> we have
					TraversalNode newBefore = new TraversalNode(node, node.node.firstKmerStartPosition(), existing.node.firstKmerStartPosition() - 1);
					memoized.add(newBefore);
					frontier.add(newBefore);
				}
				if (newStartPosition > node.node.firstKmerEndPosition()) {
					// existing node is better than us for all remaining starting position
					// -> nothing more to do
					return;
				} else {
					node = new TraversalNode(node, newStartPosition, node.node.firstKmerEndPosition());
				}
			}
		}
		memoized.add(node);
		frontier.add(node);
	}
	/**
	 * Removes returns the next node for visitation
	 * @return
	 */
	public TraversalNode pollFrontier() {
		flushInvalidFrontierHead();
		return frontier.poll();	
	}
	/**
	 * Removes returns the next node for visitation
	 * @return
	 */
	public TraversalNode peekFrontier() {
		flushInvalidFrontierHead();
		return frontier.peek();	
	}
	private void flushInvalidFrontierHead() {
		while (!frontier.isEmpty() && !isValid(frontier.peek())) {
			frontier.poll();
		}
	}
	/**
	 * Determines whether the given node is still valid, or if it has
	 * been replaced by a better path
	 * @param node node to check
	 * @return true if the given node is the highest weighted path over the given
	 * kmer position interval, false if a superior path has been found.
	 */
	private boolean isValid(TraversalNode node) {
		return memoized.contains(node);
	}
	public boolean sanityCheck() {
		// memoized scores must be unique 
		TraversalNode last = null;
		for (TraversalNode n : memoized) {
			if (last != null) {
				assert(!(last.node.kmer(0) == n.node.kmer(0) && IntervalUtil.overlapsClosed(
						last.node.firstKmerStartPosition(), last.node.firstKmerStartPosition(),
						n.node.firstKmerStartPosition(), n.node.firstKmerStartPosition())));
			}
			n = last;
		}
		return true;
	}
	private static Ordering<TraversalNode> MemoizedNodeByKmerStartEndScore = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return ComparisonChain.start()
					.compare(left.node.node().kmer(0), right.node.node().kmer(0))
					.compare(left.node.firstKmerStartPosition(), right.node.firstKmerStartPosition())
					.compare(left.node.firstKmerEndPosition(), right.node.firstKmerEndPosition())
					.compare(left.score, right.score)
					.result();
		}
	};
	private static Ordering<TraversalNode> MemoizedNodeByEnd = new Ordering<TraversalNode>() {
		@Override
		public int compare(TraversalNode left, TraversalNode right) {
			return Ints.compare(left.node.firstKmerEndPosition() + left.node.length(), right.node.firstKmerEndPosition() + right.node.length());
		}
	};
}
