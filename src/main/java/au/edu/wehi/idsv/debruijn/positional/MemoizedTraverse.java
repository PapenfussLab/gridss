package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import au.edu.wehi.idsv.util.IntervalUtil;
import it.unimi.dsi.fastutil.ints.AbstractInt2ObjectSortedMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;

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
	private final IdentityHashMap<KmerPathNode, AbstractInt2ObjectSortedMap<TraversalNode>> memoized = new IdentityHashMap<>();
	private final SortedSet<TraversalNode> frontier = new TreeSet<>(TraversalNode.ByLastEndKmer);
	/**
	 * Removes the given node from the graph
	 * @param node
	 * @return all paths unmemoized as a result of removing this node
	 */
	public void remove(KmerPathNode node) {
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
		if (cache == null) return;
		if (cache != null) {
			for (TraversalNode tn : cache.values()) {
				onMemoizeRemove(tn);
				if (frontier.remove(tn)) {
					onFrontierRemove(tn);
				}
			}
		}
		for (KmerPathNode next : node.next()) {
			// need to remove based on parent KmerPathNode, not parent
			// TraversalNode due to reference start/end edge case:
			// reference nodes can start a new path from a TraversalNode
			// which is not memoized (since a reference node can both
			// start and end separate, unconnected paths).
			unmemoize(node, next);
		}
	}
	/**
	 * Unmemoizes all paths involving the given node
	 * @param node  node to removed memoized paths from
	 */
	private void unmemoize(TraversalNode node) { 
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
		if (cache == null) return;
		TraversalNode removedNode = cache.remove(node);
		assert(removedNode == node);
		onMemoizeRemove(node);
		if (frontier.remove(node)) {
			onFrontierRemove(node);
		}
		for (KmerPathNode next : node.node.node().next()) {
			unmemoize(node, next);
		}
	}
	/**
	 * Unmemoizes all paths at the given node with the given parent
	 * @param parent parent of paths to remove
	 * @param node node to removed memoized paths from
	 * @param removedList nodes removed from memoization
	 */
	private void unmemoize(TraversalNode parent, KmerPathNode node) {
		// TODO Stack<> based version of this function so we don't overflow our actual stack
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
		if (cache == null) return;
		boolean found = true;
		// can't use a simple for loop since we may come back to ourself which will invalidate the iterator
		while (found) {
			found = false;
			// skip cached values that end before we start
			for (TraversalNode child : cache.tailMap(parent.node.firstStart()).values()) {
				if (child.parent == parent) {
					unmemoize(child);
					found = true;
					break;
				}
			}
		}
	}
	/**
	 * Unmemoizes all paths at the given node with the given parent
	 * @param parent parent of paths to remove
	 * @param node node to removed memoized paths from
	 * @param removedList nodes removed from memoization
	 */
	private void unmemoize(KmerPathNode parent, KmerPathNode node) {
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
		if (cache == null) return;
		boolean found = true;
		while (found) {
			found = false;
			for (TraversalNode child : cache.tailMap(parent.firstStart()).values()) {
				if (child.parent.node.node() == parent) { // only line not identical to function above
					unmemoize(child);
					found = true;
					break;
				}
			}
		}
	}
	/**
	 * Memoized paths for the given node
	 * @param node node
	 * @return Memoized best paths
	 */
	public Collection<TraversalNode> memoized(KmerPathNode node) {
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
		if (cache == null) return Collections.emptyList();
		return cache.values();
	}
	/**
	 * Memoize the given score for the given position
	 * @param score initial score
	 * @param start starting node
	 * @return previously memoized values invalidated by this memoization  
	 */
	public void memoize(TraversalNode node) {
		KmerPathSubnode sn = node.node;
		KmerPathNode pn = sn.node();
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(pn);
		if (cache == null) {
			cache = new Int2ObjectRBTreeMap<TraversalNode>();
			memoized.put(pn, cache);
		}
		memoize_treemap(node, cache);
	}
	private void memoize_treemap(TraversalNode node, AbstractInt2ObjectSortedMap<TraversalNode> cache) {
		KmerPathSubnode sn = node.node;
		// skip cached values that end before we start
		Iterator<TraversalNode> it = cache.tailMap(sn.firstStart()).values().iterator();
		List<TraversalNode> addlist = null;
		while (it.hasNext()) {
			TraversalNode existing = it.next();
			KmerPathSubnode existingsn = existing.node;
			assert(existingsn.firstEnd() >= sn.firstStart()); // should have been skipped in the initial lookup
			if (existingsn.firstStart() > sn.firstEnd()) break;
			
			// ok, so now we know the nodes overlap
			assert(existingsn.firstKmer() == sn.firstKmer() && IntervalUtil.overlapsClosed(existingsn.firstStart(), existingsn.firstEnd(), sn.firstStart(), sn.firstEnd()));
			if (node.score > existing.score) {
				// remove existing node in overlapping interval
				it.remove();
				boolean inFrontier = frontier.remove(existing);
				onMemoizeRemove(existing);
				if (inFrontier) {
					onFrontierRemove(existing);
				}
				
				if (existingsn.firstStart() < node.node.firstStart()) {
					// still valid in earlier interval
					TraversalNode existingBefore = new TraversalNode(existing, existingsn.firstStart(), node.node.firstStart() - 1);
					if (addlist == null) addlist = new ArrayList<TraversalNode>(4);
					addlist.add(existingBefore);
					if (inFrontier) {
						frontier.add(existingBefore);
						onFrontierAdd(existingBefore);
					}
				}
				if (existingsn.firstEnd() > node.node.firstEnd()) {
					// still valid in later interval
					TraversalNode existingAfter = new TraversalNode(existing, node.node.firstEnd() + 1, existingsn.firstEnd());
					if (addlist == null) addlist = new ArrayList<TraversalNode>(4);
					addlist.add(existingAfter);
					if (inFrontier) {
						frontier.add(existingAfter);
						onFrontierAdd(existingAfter);
					}
				}
			} else { // we're not better than existing node
				int newStartPosition = existingsn.firstEnd() + 1;
				if (node.node.firstStart() < existingsn.firstStart()) {
					// start before this node -> we are the best
					// path for the interval before the given position
					TraversalNode newBefore = new TraversalNode(node, node.node.firstStart(), existingsn.firstStart() - 1);
					if (addlist == null) addlist = new ArrayList<TraversalNode>(4);
					addlist.add(newBefore);
					frontier.add(newBefore);
					onFrontierAdd(newBefore);
				}
				if (newStartPosition > node.node.firstEnd()) {
					// existing node is better than us for all remaining starting position
					// -> don't memoize node
					node = null;
					break;
				} else {
					node = new TraversalNode(node, newStartPosition, node.node.firstEnd());
					sn = node.node;
				}
			}
		}
		if (node != null) {
			add(cache, node);
		}
		if (addlist != null) {
			for (TraversalNode n : addlist) {
				add(cache, n);
			}
		}
	}
	private void add(AbstractInt2ObjectSortedMap<TraversalNode> cache, TraversalNode node) {
		cache.put(node.node.firstEnd(), node);
		frontier.add(node);
		onMemoizeAdd(node);
	}
	/**
	 * Called when a TraversalNode is removed from memoization,
	 * including when an alternate higher scoring path is found
	 * over part of the interval
	 * @param tn
	 */
	protected void onMemoizeRemove(TraversalNode tn) {
	}
	/**
	 * Called when a TraversalNode is added as an optimal path
	 * @param tn
	 */
	protected void onMemoizeAdd(TraversalNode tn) {
	}
	protected void onFrontierAdd(TraversalNode tn) {
	}
	protected void onFrontierRemove(TraversalNode tn) {
	}
	/**
	 * Removes returns the next node for visitation
	 * @return
	 */
	public TraversalNode pollFrontier() {
		TraversalNode head = frontier.first();
		frontier.remove(head);
		onFrontierRemove(head);
		return head;
	}
	/**
	 * Removes returns the next node for visitation
	 * @return
	 */
	public TraversalNode peekFrontier() {
		return frontier.first();	
	}
	/**
	 * Chcks if the frontier is empty
	 * @return
	 */
	public boolean isEmptyFrontier() {
		return frontier.isEmpty();
	}
	/**
	 * Adds the given node back into the frontier for revisitation
	 * @param node node to recalculate
	 */
	public void addFrontier(TraversalNode node) {
		assert(memoized.containsKey(node.node.node()));
		assert(memoized.get(node.node.node()).containsValue(node));
		frontier.add(node);
		onFrontierAdd(node);
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
