package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;

import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

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
	private final IdentityHashMap<KmerPathNode, Int2ObjectRBTreeMap<MemoizedTraversalNode>> memoized = new IdentityHashMap<KmerPathNode, Int2ObjectRBTreeMap<MemoizedTraversalNode>>();
	private final PriorityQueue<MemoizedTraversalNode> frontier = new PriorityQueue<MemoizedTraversalNode>(MemoizedTraversalNode.ByLastEnd);
	/**
	 * Memoize the given score for the given position
	 * @param score initial score
	 * @param start starting node
	 */
	public void memoize(MemoizedTraversalNode node) {
		KmerPathSubnode sn = node.node;
		KmerPathNode pn = sn.node();
		Int2ObjectRBTreeMap<MemoizedTraversalNode> cache = memoized.get(pn);
		if (cache == null) {
			cache = new Int2ObjectRBTreeMap<MemoizedTraversalNode>();
			memoized.put(pn, cache);
		}
		memoize_treemap(node, cache);
	}
	private void memoize_treemap(MemoizedTraversalNode node, Int2ObjectRBTreeMap<MemoizedTraversalNode> cache) {
		KmerPathSubnode sn = node.node;
		// skip cached values that end before we start
		Iterator<MemoizedTraversalNode> it = cache.tailMap(sn.firstStart()).values().iterator();
		List<MemoizedTraversalNode> addlist = null;
		while (it.hasNext()) {
			MemoizedTraversalNode existing = it.next();
			KmerPathSubnode existingsn = existing.node;
			assert (existingsn.firstEnd() >= sn.firstStart()); // should have skipped in the initial lookup
			if (existingsn.firstStart() > sn.firstEnd()) break;
			
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
					if (addlist == null) addlist = new ArrayList<MemoizedTraversalNode>(4);
					addlist.add(existingBefore);
					if (inFrontier) frontier.add(existingBefore);
				}
				if (existingsn.firstEnd() > node.node.firstEnd()) {
					// still valid in earlier interval
					MemoizedTraversalNode existingAfter = new MemoizedTraversalNode(existing, node.node.firstEnd() + 1, existingsn.firstEnd());
					if (addlist == null) addlist = new ArrayList<MemoizedTraversalNode>(4);
					addlist.add(existingAfter);
					if (inFrontier) frontier.add(existingAfter);
				}
			} else { // we're not better than existing node
				int newStartPosition = existingsn.firstEnd() + 1;
				if (node.node.firstStart() < existingsn.firstStart()) {
					// start before this node -> we have
					MemoizedTraversalNode newBefore = new MemoizedTraversalNode(node, node.node.firstStart(), existingsn.firstStart() - 1);
					if (addlist == null) addlist = new ArrayList<MemoizedTraversalNode>(4);
					addlist.add(newBefore);
					frontier.add(newBefore);
				}
				if (newStartPosition > node.node.firstEnd()) {
					// existing node is better than us for all remaining starting position
					// -> don't memoize node
					node = null;
					break;
				} else {
					node = new MemoizedTraversalNode(node, newStartPosition, node.node.firstEnd());
					sn = node.node;
				}
			}
			
		}
		if (node != null) {
			add(cache, node);
		}
		if (addlist != null) {
			for (MemoizedTraversalNode n : addlist) {
				add(cache, n);
			}
		}
	}
	private void add(Int2ObjectRBTreeMap<MemoizedTraversalNode> cache, MemoizedTraversalNode node) {
		cache.put(node.node.firstEnd(), node);
		frontier.add(node);
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
	public int tracking_memoizedNodeCount() {
		return memoized.size();
	}
	public int tracking_frontierSize() {
		return frontier.size();
	}
}
