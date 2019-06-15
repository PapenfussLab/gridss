package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.TraversalNodeByLastEndKmerSortedSet;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker.MemoizationStats;
import com.google.common.collect.ImmutableList;
import com.google.common.io.Files;
import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.ints.AbstractInt2ObjectSortedMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectSortedMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenCustomHashSet;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;
import java.util.stream.Stream;

/**
 * Helper class to track memoization of nodes during positional graph traversal
 * @author Daniel Cameron
 *
 */
public class MemoizedTraverse {
	private static final Log log = Log.getInstance(MemoizedTraverse.class);
	/**
	 * Since a positional de Bruijn graph is a directed acyclic graph,
	 * we can calculate maximal weighted paths by a positional traverse
	 * (BFS in position space) of the graph, caching the best predecessor
	 * of each node.
	 *
	 * Int2ObjectSortedMap is sorted by TraversalNode.firstEnd
	 */
	private final IdentityHashMap<KmerPathNode, AbstractInt2ObjectSortedMap<TraversalNode>> memoized = new IdentityHashMap<>();
	// TODO: track anchored and unanchored paths in different frontiers - only call unanchored when no anchored paths nearby
	private final SortedSet<TraversalNode> frontier = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new TraversalNodeByLastEndKmerSortedSet(16) : new TreeSet<>(TraversalNode.ByLastEndKmer);
	private final MemoizationStats stats = new MemoizationStats();
	/**
	 * Removes all given nodes from the graph
	 * @param nodes nodes to remove
	 */
	public void remove(Set<KmerPathNode> nodes) {
		// bulk remove nodes in removal set
		int initialSize = memoized.size();
		Collection<TraversalNode> tns = new ArrayList<>(2 * nodes.size());
		Set<KmerPathNode> children = new ObjectOpenCustomHashSet<KmerPathNode>(new KmerNodeUtil.HashByLastEndKmer<KmerPathNode>());
		for (KmerPathNode node : nodes) {
			if (node == null) {
				if (!MessageThrottler.Current.shouldSupress(log, "removal of null KmerPathNode")) {
					log.error("Sanity check failure: removal of KmerPathNode (null)");
				}
				continue;
			}
			AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.remove(node);
			if (cache == null) {
				if (!MessageThrottler.Current.shouldSupress(log, "removal of unmemoized nodes")) {
					log.error(String.format("Sanity check failure: %s not memoized", node));
				}
			} else if (!cache.isEmpty()) {
				tns.addAll(cache.values());
				children.addAll(node.next());
			}
		}
		children.removeAll(nodes);
		onMemoizeRemove(tns);
		frontier.removeAll(tns);
		onFrontierRemove(tns);
		
		// bulk remove child paths
		Collection<TraversalNode> childPaths = removeChildPaths(children, nodes);
		int descendentCount = childPaths.size();
		onMemoizeRemove(childPaths);
		frontier.removeAll(childPaths);
		onFrontierRemove(childPaths);
		
		// Individually remove any descendant of the children
		// (bulk KmerPathNode-based removal code requires better data structure)
		int frontierResetCount = 0;
		Stack<TraversalNode> callStack = new Stack<>();
		for (TraversalNode childtn : childPaths) {
			frontierResetCount += unmemoize(childtn, callStack, true);
		}
		while (!callStack.isEmpty()) {
			frontierResetCount += unmemoize(callStack.pop(), callStack, false);
			descendentCount++;
		}
		
		stats.nodes = initialSize;
		stats.removed = nodes.size();
		stats.pathsRemoved = tns.size();
		stats.descendentPathsRemoved = descendentCount;
		stats.pathsReset = frontierResetCount;
		if (Defaults.SANITY_CHECK_MEMOIZATION) {
			assert(sanityCheckAreRemoved(nodes));
			assert(sanityCheck());
		}
	}
	/**
	 * Finds all child paths coming from any of the given parents
	 * @param toCheck nodes to check for children 
	 * @param parents parents to identify
	 */
	private Collection<TraversalNode> removeChildPaths(Iterable<KmerPathNode> toCheck, Set<KmerPathNode> parents) {
		Collection<TraversalNode> matches = new ArrayList<>();
		for (KmerPathNode node : toCheck) {
			AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
			if (cache == null || cache.isEmpty()) continue;
			Iterator<TraversalNode> it = cache.values().iterator();
			while (it.hasNext()) {
				TraversalNode tn = it.next();
				if (tn.parent != null && parents.contains(tn.parent.node.node())) {
					it.remove();
					matches.add(tn);
				}
			}
		}
		return matches;
	}
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
		Stack<TraversalNode> callStack = new Stack<TraversalNode>();
		callStack.addAll(cache.values());
		while (!callStack.isEmpty()) {
			unmemoize(callStack.pop(), callStack, false);
		}
		assert(cache.size() == 0);
		memoized.remove(node);
		if (Defaults.SANITY_CHECK_MEMOIZATION) {
			assert(sanityCheckAreRemoved(ImmutableList.of(node)));
			assert(sanityCheck());
		}
	}
	/**
	 * Unmemoizes the given path
	 * @param tn path to remove
	 * @param callStack explicit call stack to prevent stack overflow and
	 * iterator invalidation issues caused by earlier recursive implementation
	 * @param alreadyRemoved indicating whether the node has already been removed from
	 * the memoization and frontier data structures
	 * @return frontier reset count
	 */
	private int unmemoize(TraversalNode tn, Stack<TraversalNode> callStack, boolean alreadyRemoved) {
		if (!alreadyRemoved) {
			if (memoized.get(tn.node.node()).remove(tn.node.firstEnd()) == null) {
				// already processed this TraversalNode
				return 0;
			}
			onMemoizeRemove(tn);
			if (frontier.remove(tn)) {
				onFrontierRemove(tn);
			}
		}
		int frontierResetCount = addAlternatePathsToFrontier(tn);
		// check if this path continues on to any children
		for (KmerPathNode child : tn.node.node().next()) {
			AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(child);
			if (cache != null) {
				// skip values that end before we start
				for (TraversalNode childtn : cache.tailMap(tn.node.lastStart() + 1).values()) {
					// can't use reference equality since
					// the parent node could have been split
					// on an unrelated path.
					// Also, reference nodes nodes are special cases
					// and are memoized only as terminal nodes
					// with a starting node traversal object recreated
					// for each traversal
					// need to compare against the child internal since childtn.parent
					// could have an outdated overlapping interval whilest the actual
					// memoized parent interval was split and does not require updating
					if (childtn.node.firstStart() > tn.node.lastEnd() + 1) {
						break;
					}
					if (childtn.parent != null && childtn.parent.node.node() == tn.node.node()) {
						callStack.push(childtn);
					}
				}
			}
		}
		return frontierResetCount;
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
	private int addAlternatePathsToFrontier(TraversalNode tn) {
		int count = 0;
		KmerPathNode parent = tn.parent == null ? null : tn.parent.node.node();
		for (KmerPathNode prev : tn.node.node().prev()) {
			if (prev != parent) {
				int parentLength = prev.length();
				AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(prev);
				if (cache != null) {
					for (TraversalNode altParent : cache.tailMap(tn.node.firstStart() - parentLength).values()) {
						if (altParent.node.lastStart() + 1 > tn.node.firstEnd()) {
							break;
						}
						addFrontier(altParent);
					}
				}
			}
		}
		return count;
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
	 * @param node starting node
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
		assert(!cache.isEmpty());
		if (Defaults.SANITY_CHECK_MEMOIZATION) {
			//sanityCheck();
		}
	}
	private void memoize_fastutil_sortedmap(TraversalNode node, AbstractInt2ObjectSortedMap<TraversalNode> cache) {
		KmerPathSubnode sn = node.node;
		// skip cached values that end before we start
		Iterator<TraversalNode> it = cache.tailMap(sn.firstStart()).values().iterator();
		List<TraversalNode> addlist = null; // need to delay adding to cache until after our iterator is complete (so we don't invalidate the iterator midway)
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
			cache.put(node.node.firstEnd(), node);
			onMemoizeAdd(node);
			frontier.add(node);
			onFrontierAdd(node);
		}
		if (addlist != null) {
			for (TraversalNode n : addlist) {
				cache.put(n.node.firstEnd(), n);
				onMemoizeAdd(n);
			}
		}
	}
	/**
	 * Called when a TraversalNode is removed from memoization,
	 * including when an alternate higher scoring path is found
	 * over part of the interval
	 * @param tn
	 */
	protected void onMemoizeRemove(TraversalNode tn) {
	}
	protected void onMemoizeRemove(Collection<TraversalNode> tns) {
		for (TraversalNode tn : tns) {
			onMemoizeRemove(tn);
		}
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
	protected void onFrontierRemove(Collection<TraversalNode> tns) {
		for (TraversalNode tn : tns) {
			onFrontierRemove(tn);
		}
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
	 * Returns the memoized predecessor of the given path
	 * @return predecessor node, null if no single predecessor is defined for the entire interval 
	 */
	public KmerPathNode getParent(KmerPathNode node, int start, int end) {
		AbstractInt2ObjectSortedMap<TraversalNode> cache = memoized.get(node);
		Iterator<TraversalNode> it = cache.tailMap(start).values().iterator();
		if (it.hasNext()) {
			TraversalNode existing = it.next();
			if (existing.node.firstStart() <= start && existing.node.firstEnd() >= end) {
				if (existing.parent != null) {
					return existing.parent.node.node();
				}
			}
		}
		return null;
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
	public MemoizationStats tracking_lastRemoval() {
		return stats; 
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
				assert(tn.sanityCheck());
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
	public boolean sanityCheckAreRemoved(Collection<KmerPathNode> removed) {
		for (Entry<KmerPathNode, AbstractInt2ObjectSortedMap<TraversalNode>> entry : memoized.entrySet()) {
			for (int position : entry.getValue().keySet()) {
				TraversalNode tn = entry.getValue().get(position);
				sanityCheckDoesNotContain(tn, removed);
			}
		}
		for (TraversalNode tn : frontier) {
			sanityCheckDoesNotContain(tn, removed);
		}
		return true;
	}
	private void sanityCheckDoesNotContain(TraversalNode tn, Collection<KmerPathNode> removed) {
		assert(!removed.contains(tn.node.node()));
		if (tn.parent != null) {
			sanityCheckDoesNotContain(tn.parent, removed);
		}
	}
}
