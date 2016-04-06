package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.AbstractInt2ObjectSortedMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeSet;
import java.util.stream.Stream;

import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.io.Files;

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
	 * 
	 * Memoized paths involving the removed node must be removed,
	 * possibly resulting in a new starting path node.
	 * 
	 * As path weights are not known, recreation of the potential
	 * starting nodes of the immediate children is left to the
	 * caller
	 * 
	 * @param node node to remove
	 */
	public void remove(KmerPathNode node) {
		assert(node.isValid());
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
		if (cache == null) return;
		// explicit stack to prevent actual stack overflow
		Stack<TraversalNode> callStack = new Stack<TraversalNode>();
		callStack.addAll(cache.values());
		while (!callStack.isEmpty()) {
			unmemoize(callStack.pop(), callStack);
		}
		assert(cache.size() == 0);
		memoized.remove(node);
		sanityCheck();
	}
	/**
	 * Unmemoizes the given path
	 * @param tn path to remove
	 * @param callStack explicit call stack to prevent stack overflow and
	 * iterator invalidation issues caused by earlier recursive implementation
	 */
	private void unmemoize(TraversalNode tn, Stack<TraversalNode> callStack) {
		memoized.get(tn.node.node()).remove(tn.node.firstEnd());
		onMemoizeRemove(tn);
		if (frontier.remove(tn)) {
			onFrontierRemove(tn);
		}
		addAlternatePathsToFrontier(tn);
		// check if this path continues on to any children
		for (KmerPathNode child : tn.node.node().next()) {
			for (TraversalNode childtn : memoized(child)) {
				// unfortunately the expect check of
				// if (childtn.parent == tn) {
				// doesn't work since reference nodes
				// act as separate start and ends so
				// if tn is a reference node, then the
				// parent of the child tn will not match
				// the memoized path for the ref as
				// that path refers to the earlier path
				// terminating at the ref
				if (childtn.parent != null && childtn.parent.node.node() == tn.node.node() &&
						IntervalUtil.overlapsClosed(childtn.parent.node.firstStart(), childtn.parent.node.firstEnd(),
								tn.node.firstStart(), tn.node.firstEnd())
						) {
					callStack.push(childtn);
				}
			}
		}
	}
	/**
	 * Adds alternate paths to the given memoized path to
	 * the frontier.
	 * 
	 * When a memoized path is removed, the best path over the
	 * interval in which that path was the best must be recalculated.
	 * This can be done by adding all alternate paths overlapping
	 * the removed path to the frontier.
	 * 
	 * @param tn
	 */
	private void addAlternatePathsToFrontier(TraversalNode tn) {
		KmerPathNode parent = tn.parent == null ? null : tn.parent.node.node();
		for (KmerPathNode prev : tn.node.node().prev()) {
			if (prev != parent) {
				for (TraversalNode altParent : memoized(prev)) {
					// only recalculate if the interval for the alternate path
					// overlaps us
					if (IntervalUtil.overlapsClosed(tn.node.firstStart(), tn.node.firstEnd(),
							altParent.node.lastStart() + 1, altParent.node.lastEnd() + 1)) {
						addFrontier(altParent);
					}
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
		memoize_fastutil_sortedmap(node, cache);
		sanityCheck();
	}
	private void memoize_fastutil_sortedmap(TraversalNode node, AbstractInt2ObjectSortedMap<TraversalNode> cache) {
		KmerPathSubnode sn = node.node;
		// skip cached values that end before we start
		Iterator<TraversalNode> it = cache.tailMap(sn.firstStart()).values().iterator();
		List<TraversalNode> addlist = null;
		while (it.hasNext()) {
			TraversalNode existing = it.next();
			KmerPathSubnode existingsn = existing.node;
			assert(existingsn.firstEnd() >= sn.firstStart()); // should have been skipped in the initial lookup
			if (existingsn.firstStart() > sn.firstEnd()) {
				// We've already traversed all overlapping values
				// so there's no point iterating any further
				break;
			}
			// ok, so now we know the nodes overlap
			assert(existingsn.firstKmer() == sn.firstKmer() && IntervalUtil.overlapsClosed(existingsn.firstStart(), existingsn.firstEnd(), sn.firstStart(), sn.firstEnd()));
			if (node.score > existing.score) {
				// remove existing node in overlapping interval
				it.remove();
				onMemoizeRemove(existing);
				boolean inFrontier = frontier.remove(existing);
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
		onMemoizeAdd(node);
		frontier.add(node);
		onFrontierAdd(node);
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
		if (frontier.isEmpty()) return null;
		return frontier.first();	
	}
	/**
	 * Chcks if the frontier is empty
	 * @return true if the frontier is empty, false otherwise
	 */
	public boolean isEmptyFrontier() {
		return frontier.isEmpty();
	}
	/**
	 * Checks whether the given node has been memoized. A node in
	 * which all best TraversalNodes have been invalidated and require
	 * recalculation is still considered memoized.
	 * 
	 * @param node node to check
	 * @return true if at least one TraversalNode has been memoized for the given
	 * node, and the node has not been removed.
	 */
	public boolean isMemoized(KmerPathNode node) {
		return memoized.containsKey(node);
	}
	/**
	 * Adds the given node into the frontier for revisitation
	 * @param node node to recalculate
	 */
	public void addFrontier(TraversalNode node) {
		assert(memoized.containsKey(node.node.node()));
		assert(memoized.get(node.node.node()).containsValue(node));
		frontier.add(node);
		onFrontierAdd(node);
		sanityCheck();
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
	/**
	 * Exports the lookup table
	 * @param file
	 * @throws IOException 
	 */
	public void export(File file) throws IOException {
		StringBuilder sb = new StringBuilder("score,kmerlength,start,end,nodehash,parenthash,nodestart,nodeend,memoized,frontier\n");
		Stream.concat(frontier.stream(), memoized.values().stream().flatMap(m -> m.values().stream())).distinct().forEach(tn -> {
			sb.append(String.format("%d,%d,%d,%d,%x,%x,%d,%d,%b,%b\n",
				tn.score,
				tn.pathLength,
				tn.node.firstStart(),
				tn.node.firstEnd(),
				System.identityHashCode(tn.node.node()),
				tn.parent == null ? 0 : System.identityHashCode(tn.parent.node.node()),
				tn.node.node().firstStart(),
				tn.node.node().firstEnd(),
				memoized.containsKey(tn.node.node()) && memoized.get(tn.node.node()).containsValue(tn),
				frontier.contains(tn)));
		});
		Files.write(sb.toString().getBytes(), file);
	}
	public boolean sanityCheck() {
		for (Entry<KmerPathNode, AbstractInt2ObjectSortedMap<TraversalNode>> entry : memoized.entrySet()) {
			KmerPathNode node = entry.getKey();
			assert(node.isValid());
			for (int position : entry.getValue().keySet()) {
				assert(position >= node.firstStart());
				assert(position <= node.firstEnd());
				TraversalNode tn = entry.getValue().get(position);
				assert(tn.node.node() == node);
				assert(tn.node.firstEnd() == position);
			}
		}
		for (TraversalNode tn : frontier) {
			assert(memoized.containsKey(tn.node.node()));
			assert(memoized.get(tn.node.node()).get(tn.node.firstEnd()) == tn);
			assert(memoized.get(tn.node.node()).containsValue(tn));
		}
		return true;
	}
}
