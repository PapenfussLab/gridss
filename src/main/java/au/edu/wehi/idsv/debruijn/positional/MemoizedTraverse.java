package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.List;
import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.TreeSet;

import au.edu.wehi.idsv.util.IntervalUtil;

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
	private final NavigableSet<TraversalNode> memoized = new TreeSet<TraversalNode>(TraversalNode.ByKmerStartEnd);
	private final PriorityQueue<TraversalNode> frontier = new PriorityQueue<TraversalNode>(1024, TraversalNode.ByEnd);
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
		List<TraversalNode> toAdd = new ArrayList<TraversalNode>(4);
		List<TraversalNode> toRemove = new ArrayList<TraversalNode>(4);
		for (TraversalNode existing : couldOverlap) {
			if (existing.node.firstKmer() > node.node.firstKmer()) break;
			if (existing.node.firstKmer() < node.node.firstKmer()) continue;
			if (existing.node.firstEnd() < node.node.firstStart()) continue;
			if (existing.node.firstStart() > node.node.firstEnd()) break;
			
			// we overlap an existing path
			assert(existing.node.firstKmer() == node.node.firstKmer() && IntervalUtil.overlapsClosed(
					existing.node.firstStart(), existing.node.firstEnd(),
					node.node.firstStart(), node.node.firstEnd()));
			
			// ok, so now we know the nodes overlap
			if (node.score > existing.score) {
				// remove existing node in overlapping interval
				toRemove.add(existing);
				if (existing.node.firstStart() < node.node.firstStart()) {
					// still valid in earlier interval
					TraversalNode existingBefore = new TraversalNode(existing, existing.node.firstStart(), node.node.firstStart() - 1);
					toAdd.add(existingBefore);
				}
				if (existing.node.firstEnd() > node.node.firstEnd()) {
					// still valid in earlier interval
					TraversalNode existingAfter = new TraversalNode(existing, node.node.firstEnd() + 1, existing.node.firstEnd());
					toAdd.add(existingAfter);
				}
			} else { // existing node scores better than us
				int newStartPosition = existing.node.firstEnd() + 1;
				if (node.node.firstStart() < existing.node.firstStart()) {
					// start before this node -> we have
					TraversalNode newBefore = new TraversalNode(node, node.node.firstStart(), existing.node.firstStart() - 1);
					toAdd.add(newBefore);
				}
				if (newStartPosition > node.node.firstEnd()) {
					// existing node is better than us for all remaining starting position
					// -> don't memoize node
					node = null;
					break;
				} else {
					node = new TraversalNode(node, newStartPosition, node.node.firstEnd());
				}
			}
		}
		if (node != null) {
			toAdd.add(node);
		}
		// update now that we have finished iterating over memoized
		memoized.removeAll(toRemove);
		frontier.addAll(toAdd);
		memoized.addAll(toAdd);
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
		return node != null && memoized.contains(node);
	}
	@Override
	public String toString() {
		return String.format("%d nodes memoized, %d in frontier", memoized.size(), frontier.size());
	}
	public boolean sanityCheck() {
		// memoized scores must be unique 
		TraversalNode last = null;
		for (TraversalNode n : memoized) {
			if (last != null) {
				assert(!(last.node.firstKmer() == n.node.firstKmer() && IntervalUtil.overlapsClosed(
						last.node.firstStart(), last.node.firstStart(),
						n.node.firstStart(), n.node.firstStart())));
			}
			n = last;
		}
		return true;
	}
	public int tracking_memoizedNodeCount() {
		return memoized.size();
	}
	public int tracking_frontierSize() {
		return frontier.size();
	}
}
