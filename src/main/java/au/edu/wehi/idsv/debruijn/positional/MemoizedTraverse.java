package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map.Entry;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.Defaults;
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
	private final IdentityHashMap<KmerPathNode, List<MemoizedTraversalNode>> memoized = new IdentityHashMap<KmerPathNode, List<MemoizedTraversalNode>>();
	private final PriorityQueue<MemoizedTraversalNode> frontier = new PriorityQueue<MemoizedTraversalNode>(MemoizedTraversalNode.ByLastEnd);
	/**
	 * Memoize the given score for the given position
	 * @param score initial score
	 * @param start starting node
	 */
	public void memoize(MemoizedTraversalNode node) {
		KmerPathSubnode sn = node.node;
		KmerPathNode pn = sn.node();
		List<MemoizedTraversalNode> list = memoized.get(pn);
		if (list == null) {
			list = new LinkedList<MemoizedTraversalNode>();
			memoized.put(pn, list);
		}
		ListIterator<MemoizedTraversalNode> it = list.listIterator();
		while (it.hasNext()) {
			MemoizedTraversalNode existing = it.next();
			KmerPathSubnode existingsn = existing.node;
			if (existingsn.firstEnd() < sn.firstStart()) continue; // existing before
			if (existingsn.firstStart() > sn.firstEnd()) {
				// node is before this element
				it.previous(); // roll back so the insertion will be before this position
				break;
			}
			// ok, so now we know the nodes overlap
			assert(existingsn.firstKmer() == sn.firstKmer() && IntervalUtil.overlapsClosed(existingsn.firstStart(), existingsn.firstEnd(), sn.firstStart(), sn.firstEnd()));
			if (node.score > existing.score) {
				// remove existing node in overlapping interval
				it.remove();
				boolean inFrontier = existing.isValid();
				existing.invalidate();
				if (existingsn.firstStart() < node.node.firstStart()) {
					// still valid in earlier interval
					MemoizedTraversalNode existingBefore = new MemoizedTraversalNode(existing, existingsn.firstStart(), node.node.firstStart() - 1);
					it.add(existingBefore);
					if (inFrontier) frontier.add(existingBefore);
				}
				if (existingsn.firstEnd() > node.node.firstEnd()) {
					// still valid in earlier interval
					MemoizedTraversalNode existingAfter = new MemoizedTraversalNode(existing, node.node.firstEnd() + 1, existingsn.firstEnd());
					it.add(existingAfter);
					it.previous();
					if (inFrontier) frontier.add(existingAfter);
				}
			} else { // we're not better than existing node
				int newStartPosition = existingsn.firstEnd() + 1;
				if (node.node.firstStart() < existingsn.firstStart()) {
					// start before this node -> we have
					MemoizedTraversalNode newBefore = new MemoizedTraversalNode(node, node.node.firstStart(), existingsn.firstStart() - 1);
					it.previous();
					it.add(newBefore);
					frontier.add(newBefore);
				}
				if (newStartPosition > node.node.firstEnd()) {
					// existing node is better than us for all remaining starting position
					// -> don't memoize node
					return;
				} else {
					node = new MemoizedTraversalNode(node, newStartPosition, node.node.firstEnd());
					sn = node.node;
				}
			}
			
		}
		it.add(node);
		frontier.add(node);
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			sanityCheck();
		}
	}
	/**
	 * Removes returns the next node for visitation
	 * @return
	 */
	public MemoizedTraversalNode pollFrontier() {
		ensureValidFrontierHead();
		return frontier.poll();	
	}
	/**
	 * Removes returns the next node for visitation
	 * @return
	 */
	public MemoizedTraversalNode peekFrontier() {
		ensureValidFrontierHead();
		return frontier.peek();	
	}
	private void ensureValidFrontierHead() {
		while (!frontier.isEmpty() && !frontier.peek().isValid()) {
			frontier.poll();
		}
	}
	@Override
	public String toString() {
		return String.format("%d nodes memoized, %d in frontier", memoized.size(), frontier.size());
	}
	public boolean sanityCheck() {
		for (Entry<KmerPathNode, List<MemoizedTraversalNode>> entry : memoized.entrySet()) {
			assert(entry.getValue().stream().allMatch(n -> n.node.node() == entry.getKey()));
			assert(MemoizedTraversalNode.ByFirstStart.isOrdered(entry.getValue()));
			ArrayList<MemoizedTraversalNode> list = new ArrayList<MemoizedTraversalNode>(entry.getValue());
			list.sort(MemoizedTraversalNode.ByLastEnd);
			assert(MemoizedTraversalNode.ByFirstStart.isOrdered(list)); // non-overlapping
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
