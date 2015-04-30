package au.edu.wehi.idsv.graph;

import htsjdk.samtools.util.Log;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.SortedSet;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Queues;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

/**
 * Compressed graph in which each node is a path with no branches 
 * @author Daniel Cameron
 *
 */
public class PathGraph<T, PN extends PathNode<T>> implements WeightedDirectedGraph<PN> {
	private static Log log = Log.getInstance(PathGraph.class);
	protected final SubgraphAssemblyAlgorithmTracker<T, PN> tracker;
	private final PathNodeFactory<T, PN> factory;
	private final WeightedDirectedGraph<T> graph;
	private final List<PN> pathList = Lists.newArrayList();
	private final List<List<PN>> pathNext = Lists.newArrayList();
	private final List<List<PN>> pathPrev = Lists.newArrayList();
	private int pathCount = 0;
	/**
	 * Sanity checking field: total weight of graph should not change
	 */
	protected int expectedWeight;
	//private static Log log = Log.getInstance(DeBruijnPathGraph.class);
	public PathGraph(WeightedDirectedGraph<T> graph, Collection<T> seeds, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		assert(graph != null);
		assert(seeds != null);
		assert(factory != null);
		this.graph = graph;
		this.factory = factory;
		this.tracker = tracker;
		generatePathGraph(seeds);
	}
	public PathGraph(WeightedDirectedGraph<T> graph, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		this(graph, graph.allNodes(), factory, tracker);
	}
	public int getMaxNodeId() { return pathList.size() - 1; } 
	public boolean sanityCheck() {
		assert(pathNext.size() == pathList.size());
		assert(pathPrev.size() == pathList.size());
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			int weight = 0;
			int validCount = 0;
			for (int i = 0; i < pathList.size(); i++) {
				PN path = pathList.get(i);
				if (path == null) {
					assert(pathNext.get(i) == null);
					assert(pathPrev.get(i) == null);
				} else {
					validCount++;
					assert(path.getNodeId() == i);
					assert(pathNext.get(i) != null);
					assert(pathPrev.get(i) != null);
					for (PN n : next(path)) {
						assert(prev(n).contains(path));
					}
					for (PN n : prev(path)) {
						assert(next(n).contains(path));
					}
					weight += path.weight();
				}
			}
			assert(pathCount == validCount);
			assert(weight == expectedWeight);
		}
		return true;
	}
	public WeightedDirectedGraph<T> getGraph() {
		return graph;
	}
	public int getPathCount() { return pathCount; }
	public Iterable<PN> getPaths() { return Iterables.filter(pathList, Predicates.notNull()); }
	/**
	 * generates a path graph for the subgraph reachable from the given seed kmer
	 * @param seed
	 */
	private void generatePathGraph(Collection<T> seeds) {
		int kmersTraversed = 0;
		pathList.clear();
		pathNext.clear();
		pathPrev.clear();
		expectedWeight = 0;
		Queue<T> frontier = new ArrayDeque<T>();
		frontier.addAll(seeds);
		
		Map<T, PN> pathStart = Maps.newHashMap();
		Set<T> visited = Sets.newHashSet();
		while (!frontier.isEmpty()) {
			T node = frontier.poll();
			if (visited.contains(node)) continue;
			PN path = factory.createPathNode(traverseBranchless(node));
			kmersTraversed += path.length();
			pathStart.put(path.first(), path);
			visited.addAll(path.getPath());
			for (T adj : graph.prev(path.first())) {
				if (!visited.contains(adj)) {
					frontier.add(adj);
				}
			}
			for (T adj : graph.next(path.last())) {
				if (!visited.contains(adj)) {
					frontier.add(adj);
				}
			}
			// Add path to graph
			addNode(path);
			expectedWeight += path.weight();
		}
		int edges = 0;
		for (PN path : pathList) {
			assert(path != null);
			// construct edges
			List<T> nextKmers = graph.next(path.last());
			for (T adj : nextKmers) {
				assert(adj != null);
				PN next = pathStart.get(adj);
				assert(next != null);
				addEdge(path, next);
				edges++;
			}
		}
		assert(sanityCheck());
		tracker.generatePathGraph(kmersTraversed, pathCount, edges);
		shrink();
	}
	/**
	 * Returns all paths that follow the given path
	 * @param path path
	 * @return successor paths
	 */
	public List<PN> next(PN path) {
		return pathNext.get(path.getNodeId());
	}
	/**
	 * Returns all paths that follow the given path
	 * @param path path
	 * @return preceeding paths
	 */
	public List<PN> prev(PN path) {
		return pathPrev.get(path.getNodeId());
	}
	/**
	 * Returns all paths that connect to the given path
	 * @param path path
	 * @return adjacent paths
	 */
	public List<PN> adjPath(PN path) {
		List<PN> result = Lists.newArrayList(next(path));
		result.addAll(prev(path));
		return result;
	}
	/**
	 * Shrinks the graph to its minimal representation by
	 * - merging adjacent paths containing no other branches 
	 * - removing self-intersecting edges
	 * @return reduction in graph size
	 */
	public int shrink() {
		int edgesRemoved = 0;
		int nodesCollapsed = 0;
		// Keep processing until we can't anymore
		boolean changedThisIteration = true;
		while (changedThisIteration) {
			changedThisIteration = false;
			for (PN path : getPaths()) {
				if (next(path).contains(path)) {
					removeEdge(path, path);
					changedThisIteration = true;
					edgesRemoved++;
					break;
				}
				if (mergeWithNext(path)) {
					changedThisIteration = true;
					edgesRemoved++;
					nodesCollapsed++;
					break;
				}
			}
		}
		if (edgesRemoved > 0) assert(sanityCheck());
		tracker.shrink(edgesRemoved, nodesCollapsed);
		return nodesCollapsed;
	}
	/**
	 * Adds the given node to the graph
	 * @param node
	 */
	public void addNode(PN node) {
		assert(node.getNodeId() < 0);
		node.setNodeId(pathList.size());
		pathList.add(node);
		pathNext.add(new ArrayList<PN>(4));
		pathPrev.add(new ArrayList<PN>(4));
		pathCount++;
	}
	/**
	 * Removes the given node from the graph
	 * @param node node to remove
	 */
	public void removeNode(PN node) {
		int nodeId = node.getNodeId();
		assert(nodeId >= 0);
		assert(pathList.get(nodeId) != null);
		assert(pathNext.get(nodeId) != null);
		assert(pathPrev.get(nodeId) != null);
		assert(pathNext.get(nodeId).size() == 0);
		assert(pathPrev.get(nodeId).size() == 0);
		pathList.set(nodeId, null);
		pathNext.set(nodeId, null);
		pathPrev.set(nodeId, null);
		node.setNodeId(-1);
		pathCount--;
	}
	/**
	 * Adds an edge between the two nodes
	 * @param from source node
	 * @param to target node
	 */
	public void addEdge(PN from, PN to) {
		assert(from != null);
		assert(to != null);
		assert(from.getNodeId() >= 0);
		assert(to.getNodeId() >= 0);
		assert(pathList.get(from.getNodeId()) == from);
		assert(pathList.get(to.getNodeId()) == to);
		assert(!next(from).contains(to));
		assert(!prev(to).contains(from));
		next(from).add(to);
		prev(to).add(from);
	}
	/**
	 * Removes the edge between the two nodes
	 * @param from source node
	 * @param to target node
	 */
	public void removeEdge(PN from, PN to) {
		assert(from.getNodeId() >= 0);
		assert(to.getNodeId() >= 0);
		assert(pathList.get(from.getNodeId()) == from);
		assert(pathList.get(to.getNodeId()) == to);
		assert(next(from).contains(to));
		assert(prev(to).contains(from));
		next(from).remove(to);
		prev(to).remove(from);
	}
	@Override
	public String toString(Iterable<? extends PN> path) {
		return graph.toString(Iterables.concat(Iterables.transform(path, new Function<PN, Iterable<T>>() {
			@Override
			public Iterable<T> apply(PN input) {
				return input.getPath();
			}
		})));
	}
	@Override
	public int getWeight(PN node) {
		return node.weight();
	}
	@Override
	public Collection<PN> allNodes() {
		return Lists.newArrayList(getPaths());
	}
	/**
	 * Replaces edges pointing to the given node with the given node
	 * @param oldTo node to remove incoming edges from
	 * @param newTo node to add incoming edges to
	 */
	public void replaceIncomingEdges(PN oldTo, PN newTo) {
		assert(oldTo != newTo);
		for (PN prev : prev(oldTo)) {
			listReplace(next(prev), oldTo, newTo);
		}
		this.listAdd(prev(newTo), prev(oldTo));
		prev(oldTo).clear();
	}
	/**
	 * Replaces edges pointing from given node with edges from given node
	 * @param oldFrom node to remove outgoing edges from
	 * @param newFrom node to add outgoing edges to
	 */
	public void replaceOutgoingEdges(PN oldFrom, PN newFrom) {
		assert(oldFrom != newFrom);
		for (PN next : next(oldFrom)) {
			listReplace(prev(next), oldFrom, newFrom);
		}
		this.listAdd(next(newFrom), next(oldFrom));
		this.next(oldFrom).clear();
	}
	/**
	 * Attempts to merge this node with its successor 
	 * @param node node to merge with successor
	 * @return true if node could be merged without loss of information, false otherwise
	 */
	private boolean mergeWithNext(PN node) {
		if (next(node).size() != 1) return false;
		PN nextNode = next(node).get(0); 
		if (nextNode == node) return false; // don't collapse on self
		if (prev(nextNode).size() != 1) return false;
		if (prev(nextNode).get(0) != node) throw new IllegalStateException("Sanity check failure: missing matching prev entry for next node");
		
		// create new node
		PN newNode = factory.concatPathNodes(ImmutableList.of(node, nextNode));
		addNode(newNode);
		// hook up incoming and output
		replaceIncomingEdges(node, newNode);
		replaceOutgoingEdges(nextNode, newNode);
		// remove edge between the two as this is no longer required
		removeEdge(node, nextNode);
		removeNode(node);
		removeNode(nextNode);
		assert(sanityCheck());
		return true;
	}
	/**
	 * Collapses the given leaf into the given path
	 * @param toCollapse node to collapse
	 * @param path path to collapse into
	 * @param anchoredAtStart if true, path shares starting kmer with leaf, otherwise path shares ending kmer with leaf
	 */
	protected void collapseLeafInto(PN toCollapse, List<PN> path, boolean anchoredAtStart) {
		// Work out which kmers go to which path node
		int leafKmersRemaining = toCollapse.length();
		Deque<Integer> lengths = Queues.newArrayDeque();
		if (anchoredAtStart) {
			lengths.addLast(0); // no starting overhang
			for(PN node : path) {
				int nodeLength = Math.min(leafKmersRemaining, node.length());
				lengths.addLast(nodeLength);
				leafKmersRemaining -= nodeLength;
			}
			lengths.addLast(leafKmersRemaining);
		} else {
			// iterate through the other direction when we're going backwards
			lengths.addFirst(0);
			for(PN node : Lists.reverse(path)) {
				int nodeLength = Math.min(leafKmersRemaining, node.length());
				lengths.addFirst(nodeLength);
				leafKmersRemaining -= nodeLength;
			}
			lengths.addFirst(leafKmersRemaining);
		}
		assert(lengths.size() == path.size() + 2); // path is padded by start and end leaf overhangs
		int startLength = lengths.pollFirst();
		assert(startLength >= 0);
		if (startLength != 0) {
			assert(!anchoredAtStart);
			// leaf overhangs the start of the path
			List<PN> splitLeaf = split(toCollapse, ImmutableList.of(startLength, toCollapse.length() - startLength));
			toCollapse = splitLeaf.get(1);
			removeEdge(splitLeaf.get(0), toCollapse);
			addEdge(splitLeaf.get(0), path.get(0));
			// TODO correct split leaf kmers so the new path has a valid kmer sequence
		}
		int endLength = lengths.pollLast();
		assert(endLength >= 0);
		if (endLength != 0) {
			assert(anchoredAtStart);
			// leaf overhangs the end of the path
			List<PN> splitLeaf = split(toCollapse, ImmutableList.of(toCollapse.length() - endLength, endLength));
			toCollapse = splitLeaf.get(0);
			removeEdge(toCollapse, splitLeaf.get(1));
			addEdge(path.get(path.size() - 1), splitLeaf.get(1));
			// TODO correct split leaf kmers so the new path has a valid kmer sequence
		}
		if (anchoredAtStart) {
			removeEdge(prev(toCollapse).get(0), toCollapse);
		} else {
			removeEdge(toCollapse, next(toCollapse).get(0));
		}
		int offset = 0;
		for (int i = 0; i < path.size(); i++) {
			PN node = path.get(i);
			int nodeLength = lengths.pollFirst();
			assert(nodeLength > 0);
			node.merge(
					ImmutableList.of(toCollapse),
					offset,
					nodeLength,
					!anchoredAtStart ? node.length() - nodeLength : 0,
					graph);
			offset += nodeLength;
		}
		assert(lengths.isEmpty());
		removeNode(toCollapse);
	}
	/**
	 * The two paths share a common successor 
	 * @param pathA
	 * @param pathB
	 * @return true if a common successor exists, false otherwise
	 */
	protected boolean shareNextPaths(PN pathA, PN pathB) {
		// Can't just check if the end up at the same kmer as one of the
		// paths may have been collasped into another kmer path 
		// return KmerEncodingHelper.nextState(graph.getK(), pathA.getLast(), (byte)'A') == KmerEncodingHelper.nextState(graph.getK(), pathA.getLast(), (byte)'A');
		// Check the next states are the same
		// Since next states are ordered, we can just compare that the two lists are the same 
		return Iterables.elementsEqual(next(pathA), next(pathB));
	}
	private static int inconsistentMergePathRemainingCalls = 8;
	/**
	 * Merges the given paths together
	 * @param pathA first path to merge 
	 * @param pathB second path to merge
	 * @return true if a merge could be performed, false otherwise
	 */
	public boolean mergePaths(Iterable<PN> pathA, Iterable<PN> pathB) {
		PN problematicNode = getNodeAtDifferentPosition(pathA, pathB); 
		if (problematicNode != null) {
			if (log != null && inconsistentMergePathRemainingCalls > 0) {
				inconsistentMergePathRemainingCalls--;
				log.debug(String.format("Near %s: found similar but inconsistent paths \"%s\" and \"%s\". Both contain \"%s\"",
					problematicNode.toString(getGraph()),
					toString(pathA),
					toString(pathB),
					toString(ImmutableList.of(problematicNode))));
			}
			return false;
		}
		int weightA = PathNode.totalWeight(pathA);
		int weightB = PathNode.totalWeight(pathB);
		
		Iterable<PN> imainPath = weightA >= weightB ? pathA : pathB;
		Iterable<PN> ialtPath = weightA >= weightB ? pathB : pathA;
		
		SortedSet<Integer> breaksAt = Sets.newTreeSet();
		breaksAt.addAll(getBreaks(pathA));
		breaksAt.addAll(getBreaks(pathB));
		List<PN> mainPath = splitPathToEnsureBreaksAt(imainPath, breaksAt);
		List<PN> altPath = splitPathToEnsureBreaksAt(ialtPath, breaksAt);
		assert(mainPath.size() == altPath.size());
		for (int i = 0; i < mainPath.size(); i++) {
			mergePath(altPath.get(i), mainPath.get(i));
		}
		// The collapse could have resulted in paths
		// which are now branchless but contain
		// multiple nodes: merge those nodes together
		// (eg: A-B1 in the example collapsePaths example)
		shrink();
		assert(sanityCheck());
		return true;
	}
	/**
	 * Checks for a common node that cannot be collapsed as it occurs
	 * on both paths at different positions
	 * 
	 * @param pathA
	 * @param pathB
	 * @return first inconsistent node, null if all nodes are consistent
	 */
	protected PN getNodeAtDifferentPosition(Iterable<PN> pathA, Iterable<PN> pathB) {
		Map<PN, Integer> offsetMap = Maps.newHashMap();
		int offset = 0;
		for (PN n : pathA) {
			offsetMap.put(n, offset);
			offset += n.length();
		}
		int bOffset = 0;
		for (PN b : pathB) {
			if (offsetMap.containsKey(b)) {
				int aOffset = offsetMap.get(b);
				if (aOffset != bOffset) return b;
			}
			bOffset += b.length();
		}
		return null;
	}
	private SortedSet<Integer> getBreaks(Iterable<PN> path) {
		SortedSet<Integer> breaksAt = Sets.newTreeSet();
		int offset = 0;
		breaksAt.add(0);
		for (PN n : path) {
			offset += n.length();
			breaksAt.add(offset);
		}
		return breaksAt;
	}
	/**
	 * Splits nodes on the given path to ensure breaks exist at
	 * all the given kmer offsets
	 * @param path kmer path to break
	 * @param breaksAt kmer offset of breaks required positions 
	 * @return path including breaks
	 */
	protected List<PN> splitPathToEnsureBreaksAt(Iterable<PN> path, SortedSet<Integer> breaksAt) {
		List<PN> result = Lists.newArrayList();
		int offset = 0;
		for (PN n : path) {
			SortedSet<Integer> breaksOffsets = breaksAt.subSet(offset + 1, offset + n.length());
			List<Integer> splitLength = Lists.newArrayListWithCapacity(breaksOffsets.size() + 1);
			int lastBreak = offset;
			for (int breakPos : breaksOffsets) {
				splitLength.add(breakPos - lastBreak);
				lastBreak = breakPos;
			}
			splitLength.add(n.length() - (lastBreak - offset));
			result.addAll(split(n, splitLength));
			offset += n.length();
		}
		return result;
	}
	/**
	 * Merges the given paths together
	 * @param merge path to merge
	 * @param into path to merge into
	 */
	private void mergePath(PN merge, PN into) {
		if (merge == into) return; // nothing to do
		assert(merge.length() == into.length());
		replaceIncomingEdges(merge, into);
		replaceOutgoingEdges(merge, into);
		into.merge(ImmutableList.of(merge), getGraph());
		removeNode(merge);
	}
	/**
	 * Splits a pathnode at the given offset
	 * @param node
	 * @param lengths kmer lengths of nodes to split at
	 */
	public List<PN> split(PN node, List<Integer> lengths) {
		if (lengths.size() == 0) throw new IllegalArgumentException("Split must have at least one length");
		if (lengths.size() == 1) {
			assert(lengths.get(0) == node.length());
			return ImmutableList.of(node);
		}
		int sumlength = 0;
		for (int i = 0; i < lengths.size(); i++) {
			if (lengths.get(i) <= 0) throw new IllegalArgumentException("path lengths must be greater than zero");
			sumlength += lengths.get(i);
		}
		if (sumlength != node.length()) throw new IllegalArgumentException("path lengths must sum to node path length");
		// replace the source node with the split ones in the path graph
		List<PN> result = Lists.newArrayList();
		int offset = 0;
		for (int i = 0; i < lengths.size(); i++) {
			result.add(factory.splitPathNode(node, offset, lengths.get(i)));
			offset += lengths.get(i);
		}
		assert(offset == node.length());
		// now replace the node in the graph with the new nodes
		for (PN newNode : result) {
			addNode(newNode);
		}
		replaceIncomingEdges(node, result.get(0));
		for (int i = 0; i < result.size() - 1; i++) {
			next(result.get(i)).add(result.get(i + 1));
			prev(result.get(i + 1)).add(result.get(i));
		}
		replaceOutgoingEdges(node, result.get(result.size() - 1));
		removeNode(node);
		return result;
	}
	/**
	 * Adds the given PathNodes to the given list
	 * @param toAddTo list to add to
	 * @param toAdd list to add
	 */
	private void listAdd(List<PN> toAddTo, List<PN> toAdd) {
		for (PN node : toAdd) {
			if (!toAddTo.contains(node)) {
				toAddTo.add(node);
			}
		}
	}
	private void listReplace(List<PN> list, PN toReplace, PN replaceWith) {
		if (list.contains(toReplace)) {
			if (list.contains(replaceWith)) {
				// contains both: just remove the element to be replaced as the
				// target element is already in the list
				list.remove(toReplace);
			} else {
				// only contains toReplace, just overwrite with the new value
				list.set(list.indexOf(toReplace), replaceWith);
			}
		}
	}
	/**
	 * A bubble is a graph path diverging from a reference path
	 * (usually due to a sequencing error), then converges back to the reference
	 * without any other branches
	 * 
	 * @param path path to test
	 * @return true if the path could be a bubble, false otherwise
	 */
	public boolean isBubble(Iterable<PN> path) {
		return Iterables.all(path, new Predicate<PN>() {
			@Override
			public boolean apply(PN node) {
				return prev(node).size() == 1 && next(node).size() == 1;
			}
		});
	}
	/**
	 * Traverses source graph until a branch is found
	 * @param seed starting node
	 * @return unique path
	 */
	private LinkedList<T> traverseBranchless(T seed) {
		LinkedList<T> path = new LinkedList<T>();
		Set<T> visited = Sets.newHashSet();
		path.add(seed);
		visited.add(seed);
		for(List<T> adj = graph.next(path.getLast()); adj.size() == 1 && graph.prev(adj.get(0)).size() <= 1; adj = graph.next(path.getLast())) {
			if (visited.contains(adj.get(0))) {
				// circular contig
				break;
			}
			path.addLast(adj.get(0));
			visited.add(adj.get(0));
		}
		for(List<T> adj = graph.prev(path.getFirst()); adj.size() == 1 && graph.next(adj.get(0)).size() <= 1; adj = graph.prev(path.getFirst())) {
			if (visited.contains(adj.get(0))) {
				// circular contig
				break;
			}
			path.addFirst(adj.get(0));
			visited.add(adj.get(0));
		}
		return path;
	}
	/**
	 *  Greedily traverse a path based on the given choice ordering
	 * @param startNode starting position
	 * @param forwardChoice ordering when traversing forward
	 * @param backwardChoice ordering when traversing backward 
	 * @param excluded these paths cannot be traversed
	 * @return traversal
	 */
	public LinkedList<PN> greedyTraverse(PN startNode, Comparator<PN> forwardChoice, Comparator<PN> backwardChoice, Set<PN> excluded) {
		LinkedList<PN> path = new LinkedList<PN>();
		Set<PN> visited = Sets.newHashSet();
		if (excluded != null) visited.addAll(excluded);
		path.add(startNode);
		visited.add(startNode);
		// assemble forward
		PriorityQueue<PN> nextCandidates = new PriorityQueue<PN>(4, forwardChoice);
		for (List<PN> nodeList = next(path.getLast()); ; nodeList = next(path.getLast())) {
			nextCandidates.clear();
			for (PN node : nodeList) {
				if (visited.contains(node)) continue; // no loops
				nextCandidates.add(node);
			}
			// we're done
			if (nextCandidates.isEmpty()) break;
			PN node = nextCandidates.poll();
			path.addLast(node);
			visited.add(node);
		}
		// assemble back
		PriorityQueue<PN> prevCandidates = new PriorityQueue<PN>(4, backwardChoice);
		for (List<PN> nodeList = prev(path.getFirst()); ; nodeList = prev(path.getFirst())) {
			prevCandidates.clear();
			for (PN node : nodeList) {
				if (visited.contains(node)) continue; // no loops
				prevCandidates.add(node);
			}
			// we're done
			if (prevCandidates.isEmpty()) break;
			PN node = prevCandidates.poll();
			path.addFirst(node);
			visited.add(node);
		}
		return path;
	}
	public PN getNodeContaining(T node) {
		for (PN pn : getPaths()) {
			if (pn.contains(node)) return pn;
		}
		return null;
	}
	/**
	 * Ordering by average
	 */
	public Ordering<PN> ByMeanNodeWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Doubles.compare(o1.weight() / o1.length(), o2.weight() / o2.length());
		}
	}.reverse();
	/**
	 * Ordering by total weight of each path.
	 */
	public Ordering<PN> ByPathTotalWeight = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Ints.compare(o1.weight(), o2.weight());
		}
	};
	public Ordering<PN> ByPathTotalWeightDesc = ByPathTotalWeight.reverse();
	/**
	 * Ordering by max node weight
	 */
	public Ordering<PN> ByMaxNodeWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Ints.compare(o1.maxNodeWeight(graph), o2.maxNodeWeight(graph));
		}
	}.reverse();
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("Path graph: %d node\n", pathCount));
		for (PN x : getPaths()) {
			sb.append(x.toString(getGraph()));
			sb.append(" <-{");
			for (PN s : prev(x)) sb.append(String.format("%d,", s.getNodeId()));
			sb.append("} ->{");
			for (PN s : next(x)) sb.append(String.format("%d,", s.getNodeId()));
			sb.append("}\n");
		}
		return sb.toString();
	}
}
