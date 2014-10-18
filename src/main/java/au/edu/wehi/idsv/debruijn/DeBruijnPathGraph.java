package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.util.Log;

import java.util.ArrayDeque;
import java.util.Comparator;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.SortedSet;

import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;

import com.google.common.base.Predicate;
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
 * Compressed De Bruijn graph in which each node is a kmer path with no branches 
 * @author Daniel Cameron
 *
 */
public class DeBruijnPathGraph<T extends DeBruijnNodeBase, PN extends PathNode<T>> {
	private static Log log = Log.getInstance(DeBruijnPathGraph.class);
	/**
	 * Safety limit to prevent unbounded exponential runtime
	 * when attempt to path collapse highly collected degenerate subgraphs
	 */
	private static final int COLLAPSE_PATH_MAX_TRAVERSAL = 10 * 1024 * 1024;
	protected final PathNodeFactory<T, PN> factory;
	private final DeBruijnGraphBase<T> graph;
	protected final Set<PN> paths = Sets.newHashSet();
	protected final Map<PN, List<PN>> pathNext = Maps.newHashMap();
	protected final Map<PN, List<PN>> pathPrev = Maps.newHashMap();
	private int nodeTaversalCount;
	/**
	 * Sanity checking field: total weight of graph should not change
	 */
	protected int expectedWeight;
	//private static Log log = Log.getInstance(DeBruijnPathGraph.class);
	public DeBruijnPathGraph(DeBruijnGraphBase<T> graph, long seed, PathNodeFactory<T, PN> factory) {
		this.graph = graph;
		this.factory = factory;
		generatePathGraph(seed);
	}
	public DeBruijnPathGraph(DeBruijnGraphBase<T> graph, PathNodeFactory<T, PN> factory) {
		this.graph = graph;
		this.factory = factory;
		generatePathGraph(null);
	}
	public boolean sanityCheck() {
		int weight = 0;
		assert(pathNext.size() == paths.size());
		assert(pathPrev.size() == paths.size());
		for (PN path : paths) {
			assert(pathNext.containsKey(path));
			assert(pathPrev.containsKey(path));
			for (PN n : nextPath(path)) {
				assert(prevPath(n).contains(path));
			}
			for (PN n : prevPath(path)) {
				assert(nextPath(n).contains(path));
			}
			weight += path.getWeight();
		}
		assert(weight == expectedWeight);
		return true;
	}
	public DeBruijnGraphBase<T> getGraph() {
		return graph;
	}
	public Set<PN> getPaths() { return paths; }
	/**
	 * generates a path graph for the subgraph reachable from the given seed kmer
	 * @param seed
	 */
	private void generatePathGraph(Long seed) {
		paths.clear();
		pathNext.clear();
		pathPrev.clear();
		expectedWeight = 0;
		Queue<Long> frontier = new ArrayDeque<Long>();
		if (seed == null) {
			frontier.addAll(getGraph().getAllKmers());
		} else {
			frontier.add(seed);
		}
		Map<Long, PN> pathStart = Maps.newHashMap();
		Set<Long> visited = Sets.newHashSet();
		while (!frontier.isEmpty()) {
			long kmer = frontier.poll();
			if (visited.contains(kmer)) continue;
			PN path = factory.createPathNode(traverseBranchless(kmer), getGraph());
			pathStart.put(path.getFirst(), path);
			visited.addAll(path.getPath());
			for (long adj : getGraph().prevStates(path.getFirst(), null, null)) {
				frontier.add(adj);
			}
			for (long adj : getGraph().nextStates(path.getLast(), null, null)) {
				frontier.add(adj);
			}
			// Add path to graph
			addNode(path);
			expectedWeight += path.getWeight();
		}
		for (PN path : paths) {
			// construct edges
			List<Long> nextKmers = getGraph().nextStates(path.getLast(), null, null);
			for (long adj : nextKmers) {
				PN next = pathStart.get(adj);
				addEdge(path, next);
			}
		}
		assert(sanityCheck());
		shrink();
	}
	/**
	 * Returns all paths that follow the given path
	 * @param path path
	 * @return successor paths
	 */
	public List<PN> nextPath(PN path) {
		return pathNext.get(path);
	}
	/**
	 * Returns all paths that follow the given path
	 * @param path path
	 * @return preceeding paths
	 */
	public List<PN> prevPath(PN path) {
		return pathPrev.get(path);
	}
	/**
	 * Returns all paths that connect to the given path
	 * @param path path
	 * @return adjacent paths
	 */
	public List<PN> adjPath(PN path) {
		List<PN> result = Lists.newArrayList(nextPath(path));
		result.addAll(prevPath(path));
		return result;
	}
	/**
	 * Shrinks the graph to its minimal representation by
	 * - merging adjacent paths containing no other branches 
	 * - removing self-intersecting edges
	 * @return true if the graph was update, false if no changes were required
	 */
	public boolean shrink() {
		boolean anyChange = false;
		// Keep processing until we can't anymore
		boolean changedThisIteration = true;
		while (changedThisIteration) {
			changedThisIteration = false;
			for (PN path : paths) {
				if (nextPath(path).contains(path)) {
					removeEdge(path, path);
					changedThisIteration = true;
					break;
				}
				// TODO: collapse the path with the weakest support first
				if (mergeWithNext(path)) {
					changedThisIteration = true;
					anyChange = true;
					break;
				}
			}
		}
		if (anyChange) assert(sanityCheck());
		return anyChange;
	}
	/**
	 * Adds the given node to the graph
	 * @param node
	 */
	protected void addNode(PN node) {
		assert(!paths.contains(node));
		paths.add(node);
		pathNext.put(node, Lists.<PN>newArrayListWithExpectedSize(4));
		pathPrev.put(node, Lists.<PN>newArrayListWithExpectedSize(4));
	}
	/**
	 * Removes the given node from the graph
	 * @param node node to remove
	 */
	protected void removeNode(PN node) {
		assert(paths.contains(node));
		assert(pathNext.containsKey(node));
		assert(pathPrev.containsKey(node));
		assert(pathNext.get(node).size() == 0);
		assert(pathPrev.get(node).size() == 0);
		paths.remove(node);
		pathNext.remove(node);
		pathPrev.remove(node);
	}
	/**
	 * Adds an edge between the two nodes
	 * @param from source node
	 * @param to target node
	 */
	protected void addEdge(PN from, PN to) {
		assert(paths.contains(from));
		assert(paths.contains(to));
		assert(!nextPath(from).contains(to));
		assert(!prevPath(to).contains(from));
		nextPath(from).add(to);
		prevPath(to).add(from);
	}
	/**
	 * Removes the edge between the two nodes
	 * @param from source node
	 * @param to target node
	 */
	protected void removeEdge(PN from, PN to) {
		assert(paths.contains(from));
		assert(paths.contains(to));
		assert(nextPath(from).contains(to));
		assert(prevPath(to).contains(from));
		nextPath(from).remove(to);
		prevPath(to).remove(from);
	}
	/**
	 * Replaces edges pointing to the given node with the given node
	 * @param oldTo node to remove incoming edges from
	 * @param newTo node to add incoming edges to
	 */
	public void replaceIncomingEdges(PN oldTo, PN newTo) {
		assert(oldTo != newTo);
		for (PN prev : prevPath(oldTo)) {
			listReplace(nextPath(prev), oldTo, newTo);
		}
		this.listAdd(prevPath(newTo), prevPath(oldTo));
		prevPath(oldTo).clear();
	}
	/**
	 * Replaces edges pointing from given node with edges from given node
	 * @param oldFrom node to remove outgoing edges from
	 * @param newFrom node to add outgoing edges to
	 */
	public void replaceOutgoingEdges(PN oldFrom, PN newFrom) {
		assert(oldFrom != newFrom);
		for (PN next : nextPath(oldFrom)) {
			listReplace(prevPath(next), oldFrom, newFrom);
		}
		this.listAdd(nextPath(newFrom), nextPath(oldFrom));
		this.nextPath(oldFrom).clear();
	}
	/**
	 * Attempts to merge this node with its successor 
	 * @param node node to merge with successor
	 * @return true if node could be merged without loss of information, false otherwise
	 */
	private boolean mergeWithNext(PN node) {
		if (nextPath(node).size() != 1) return false;
		PN next = nextPath(node).get(0); 
		if (next == node) return false; // don't collapse on self
		if (prevPath(next).size() != 1) return false;
		if (prevPath(next).get(0) != node) throw new IllegalStateException("Sanity check failure: missing matching prev entry for next node");
		
		// create new node
		PN newNode = factory.createPathNode(getGraph(), ImmutableList.of(node, next));
		addNode(newNode);
		// hook up incoming and output
		replaceIncomingEdges(node, newNode);
		replaceOutgoingEdges(next, newNode);
		// remove edge between the two as this is no longer required
		removeEdge(node, next);
		removeNode(node);
		removeNode(next);
		assert(sanityCheck());
		return true;
	}
	/**
	 * Collapses all paths that differ by at most maxDifference bases
	 * into the path with the highest weight.
	 * @param maxDifference maximum differences along the path. Note that by collapsing
	 * a graph with n differences, alternate paths that share a PathNode  may have up to n differences
	 * from their original due to a differences merging the paths
	 *       E      F                    E             F
	 *        \    /                      \           / 
	 *    G - H - I - J    ->              \    K    / 
	 *   /              \                   \  / \  /        
	 *  A  -  B  -  C -  D          A - B1 - B2 - C1 - C2 - D
	 *         \   /                ^^^^^               ^^^^^   <- branches paths needing merging
	 *           K
	 *           
	 *  E-H-I-J is turned in E-B1-B2-C1-C2-F resulting in a different path
	 *  
	 *  @param bubblesOnly only collapse bubbles. No other paths will be affected
	 *  @return number of paths collapsed
	 */
	public int collapseSimilarPaths(int maxDifference, boolean bubblesOnly) {
		nodeTaversalCount = 0;
		int collapseCount = 0;
		try {
			// Keep collapsing until we can't anymore
			boolean collapsed = true;
			while (collapsed) {
				collapsed = false;
				for (PN start : paths) {
					// TODO: collapse the path with the weakest support first
					if (collapsePaths(maxDifference, bubblesOnly, start)) {
						collapseCount++;
						collapsed = true;
						break;
					}
				}
			}
		} catch (AlgorithmRuntimeSafetyLimitExceededException e) {
			log.error(e);
		}
		return collapseCount;
	}
	/**
	 * Collapses low-weight leaf branches into higher weighted paths
	 * @param maxDifference maximum base pair difference
	 * @param start starting node
	 * @return number of leaves collapsed
	 * @throws AlgorithmRuntimeSafetyLimitExceededException
	 */
	public int collapseLeaves(int maxDifference) {
		nodeTaversalCount = 0;
		int collapseCount = 0;
		try {
			int collapsedThisRound;
			do {
				collapsedThisRound = 0;
				PriorityQueue<PN> ordering = new PriorityQueue<PN>(paths.size(), ByPathTotalWeightDesc.reverse());
				ordering.addAll(paths);
				// unlike collapsing paths, collapsing leaves is a local
				// change and does not invalidate any other paths
				for (PN leaf : ordering) {
					List<PN> nextPaths = nextPath(leaf);
					List<PN> prevPaths = prevPath(leaf);
					if (nextPaths.size() == 0 && prevPaths.size() == 1) {
						if (collapseLeaf(maxDifference, leaf, prevPaths.get(0), true)) {
							collapseCount++;
							collapsedThisRound++;
						}
					} else if (nextPaths.size() == 1 && prevPaths.size() == 0) {
						if (collapseLeaf(maxDifference, leaf, nextPaths.get(0), false)) {
							collapseCount++;
							collapsedThisRound++;
						}
					}
				}
				if (collapsedThisRound > 0) {
					// leaf collapse may have resulted in nodes that can now be merged
					shrink();
				}
				// the shrink process could have made a new leaf that can now be collapsed
				// so we need to try again
				// Optimisation: only reprocess shrunk nodes 
			} while (collapsedThisRound > 0);
		} catch (AlgorithmRuntimeSafetyLimitExceededException e) {
			log.error(e);
		}
		return collapseCount;
	}
	/**
	 * 
	 * @param maxDifference maximum bases different
	 * @param leaf leaf node
	 * @param anchor anchor to start collapse
	 * @param forward anchor to leaf is in the direction of kmer traversal 
	 * @return true if the leaf was collapsed, false if not similar path found
	 * @throws AlgorithmRuntimeSafetyLimitExceededException 
	 */
	private boolean collapseLeaf(int maxDifference, PN leaf, PN anchor, boolean forward) throws AlgorithmRuntimeSafetyLimitExceededException {
		// find path from anchor that is within maxDifference from leaf
		List<PN> path = findHighestWeightSimilarPath(ImmutableList.<PN>of(), Queues.<PN>newArrayDeque(), 0, 0, maxDifference, leaf, anchor, forward);
		if (path == null || path.isEmpty()) return false;
		// compare weights of shared kmers
		int pathLength = PathNode.kmerLength(path);
		int sharedLength = Math.min(leaf.length(), pathLength);
		if (PathNode.kmerTotalWeight(ImmutableList.of(leaf), forward ? 0 : leaf.length() - sharedLength, sharedLength, getGraph()) >=
			PathNode.kmerTotalWeight(path, forward ? 0 : pathLength - sharedLength, sharedLength, getGraph())) {
			return false;
		}
		collapseLeafInto(leaf, path, forward);
		assert(sanityCheck());
		return true;
	}
	/**
	 * Collapses the given leaf into the given path
	 * @param toCollapse node to collapse
	 * @param path path to collapse into
	 * @param anchoredAtStart if true, path shares starting kmer with leaf, otherwise path shares ending kmer with leaf
	 */
	private void collapseLeafInto(PN toCollapse, List<PN> path, boolean anchoredAtStart) {
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
			removeEdge(prevPath(toCollapse).get(0), toCollapse);
		} else {
			removeEdge(toCollapse, nextPath(toCollapse).get(0));
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
	 * Finds all paths similar to the leaf path and returns the highest weighted similar path. 
	 * @param bestPath best path so far
	 * @param currentPath current search path
	 * @param currentLength length of current search path
	 * @param currentDifferences current number of search path bases different
	 * @param maxDifference maximum path differences
	 * @param leaf path to compare to
	 * @param anchor anchor path that leaf and search paths share
	 * @param traverseForward direction of traversal
	 * @return similar path with highest weight
	 * @throws AlgorithmRuntimeSafetyLimitExceededException thrown when maximum number of search nodes has been exceeded
	 */
	private List<PN> findHighestWeightSimilarPath(
			List<PN> bestPath,
			Deque<PN> currentPath,
			int currentLength,
			int currentDifferences,
			
			int maxDifference,
			PN leaf,
			PN anchor,
			boolean traverseForward) throws AlgorithmRuntimeSafetyLimitExceededException {
		if (++nodeTaversalCount >= COLLAPSE_PATH_MAX_TRAVERSAL) {
			throw new AlgorithmRuntimeSafetyLimitExceededException(String.format(
					"Leaf collapse not complete after %d node traversals. Aborting path collapse whilst processing \"%s\".",
					nodeTaversalCount,
					new String(getGraph().getBaseCalls(leaf.getPath()))));
		}
		if (currentDifferences > maxDifference) {
			return bestPath;
		}
		List<PN> nextList = traverseForward ? nextPath((currentPath.isEmpty() ? anchor : currentPath.getLast())) : prevPath((currentPath.isEmpty() ? anchor : currentPath.getFirst()));
		if (currentLength >= leaf.length()) {
			assert(currentLength == PathNode.kmerLength(currentPath));
			int leafLength = leaf.length();
			int bestLength = PathNode.kmerLength(bestPath);
			int currentSharedLength = Math.min(currentLength, leafLength);
			int bestSharedLength = Math.min(bestLength, leafLength);
			if (PathNode.kmerTotalWeight(currentPath, traverseForward ? 0 : currentLength - currentSharedLength, currentSharedLength, getGraph()) >
				PathNode.kmerTotalWeight(bestPath, traverseForward ? 0 : bestLength - bestSharedLength, bestSharedLength, getGraph())) {
				// this path has a higher total weight over the kmers shared with the leaf 
				return Lists.newArrayList(currentPath);
			} else {
				return bestPath;
			}
		}
		if (nextList.size() == 0) {
			// Don't merge a leaf into shorter paths
			// Doing so messes up consensus kmer sequence - the primary sequence
			// would no longer be a valid kmer path. We also can't try to fix the
			// consensus sequence as the kmer we replace the inconsistent kmers
			// with could already be part of our graph elsewhere
			// or, even worse, part of a different subgraph in the kmer graph! 
			return bestPath;
		}
		for (PN next : nextList) {
			if (next != leaf && !currentPath.contains(next)) {
				if (traverseForward) {
					currentPath.addLast(next);
				} else {
					currentPath.addFirst(next);
				}
				bestPath = findHighestWeightSimilarPath(
						bestPath,
						currentPath,
						currentLength + next.length(),
						traverseForward ? basesDifferent(ImmutableList.of(leaf), currentPath) : reverseBasesDifferent(ImmutableList.of(leaf), currentPath),
						maxDifference, leaf, anchor, traverseForward);
				if (traverseForward) {
					currentPath.removeLast();
				} else {
					currentPath.removeFirst();
				}
			}
		}
		return bestPath;
	}
	private boolean collapsePaths(int maxDifference, boolean bubblesOnly, PN start) throws AlgorithmRuntimeSafetyLimitExceededException {
		Deque<PN> listA = Queues.newArrayDeque();
		Deque<PN> listB = Queues.newArrayDeque();
		List<PN> next = nextPath(start);
		for (int i = 0; i < next.size() - 1; i++) {
			for (int j = i + 1; j < next.size(); j++) {
				listA.add(next.get(i));
				listB.add(next.get(j));
				if (collapsePaths(maxDifference, bubblesOnly, listA, listB, next.get(i).length(), next.get(j).length())) {
					return true;
				}
				assert(listA.size() == 1);
				assert(listB.size() == 1);
				listA.removeFirst();
				listB.removeFirst();
			}
		}
		return false;
	}
	private boolean collapsePaths(
			int differencesAllowed,
			boolean bubblesOnly,
			Deque<PN> pathA,
			Deque<PN> pathB,
			int pathALength,
			int pathBLength) throws AlgorithmRuntimeSafetyLimitExceededException {
		if (++nodeTaversalCount >= COLLAPSE_PATH_MAX_TRAVERSAL) {
			throw new AlgorithmRuntimeSafetyLimitExceededException(String.format(
					"Path collapse not complete after %d node traversals. Aborting path collapse whilst processing \"%s\" and \"%s\".",
					nodeTaversalCount,
					new String(getGraph().getBaseCalls(Lists.newArrayList(PathNode.kmerIterator(pathA)))),
					new String(getGraph().getBaseCalls(Lists.newArrayList(PathNode.kmerIterator(pathB))))));
		}
		// paths have diverged too far
		// Optimisation: cache bases compared and difference so far so only 2 PNs need to be compared
		if (basesDifferent(pathA, pathB) > differencesAllowed) return false;
		
		if (shareNextPaths(pathA.getLast(), pathB.getLast())) {
			// only accept base difference, not indels
			if (pathALength != pathBLength) return false;			
			// we have a common path!
			if (bubblesOnly && (!isBubble(pathA) && !isBubble(pathB))) return false;
			boolean collapsed = mergePaths(pathA, pathB);
			return collapsed;
		}
		// short circuit if only processing bubbles
		if (bubblesOnly && (!isBubble(pathA) && !isBubble(pathB))) return false;
		
		if (pathALength <= pathBLength) {
			for (PN nextA : nextPath(pathA.getLast())) {
				if (!pathA.contains(nextA)) {
					pathA.addLast(nextA);
					// don't continue if the node we just added would result in a problematic path merge
					if (!pathB.contains(nextA) || getNodeAtDifferentPosition(pathA, pathB) == null) {
						if (collapsePaths(differencesAllowed, bubblesOnly, pathA, pathB, pathALength + nextA.length(), pathBLength)) {
							pathA.removeLast();
							return true;
						}
					}
					pathA.removeLast();
				}
			}
		} else {
			for (PN nextB : nextPath(pathB.getLast())) {
				if (!pathB.contains(nextB)) {
					pathB.addLast(nextB);
					if (!pathA.contains(nextB) || getNodeAtDifferentPosition(pathA, pathB) == null) {
						if (collapsePaths(differencesAllowed, bubblesOnly, pathA, pathB, pathALength, pathBLength + nextB.length())) {
							pathB.removeLast();
							return true;
						}
					}
					pathB.removeLast();
				}
			}
		}
		return false; 
	}
	/**
	 * The two paths share a common successor 
	 * @param pathA
	 * @param pathB
	 * @return true if a common successor exists, false otherwise
	 */
	private boolean shareNextPaths(PN pathA, PN pathB) {
		// Can't just check if the end up at the same kmer as one of the
		// paths may have been collasped into another kmer path 
		// return KmerEncodingHelper.nextState(graph.getK(), pathA.getLast(), (byte)'A') == KmerEncodingHelper.nextState(graph.getK(), pathA.getLast(), (byte)'A');
		// Check the next states are the same
		// Since next states are ordered, we can just compare that the two lists are the same 
		return Iterables.elementsEqual(nextPath(pathA), nextPath(pathB));
	}
	/**
	 * Returns the number of bases difference between the two paths when
	 * traversing the paths backward.
	 * Bases are only compared up to the length of the shortest of the
	 * paths
	 * @param pathA
	 * @param pathB
	 * @return number of bases difference between the two paths
	 */
	protected int reverseBasesDifferent(Iterable<PN> pathA, Iterable<PN> pathB) {
		List<Long> lA = Lists.newArrayList(PathNode.kmerIterator(pathA));
		List<Long> lB = Lists.newArrayList(PathNode.kmerIterator(pathB));
		// skip the initial bases of the longer path
		Iterator<Long> itA = lA.listIterator(Math.max(0, lA.size() - lB.size()));
		Iterator<Long> itB = lB.listIterator(Math.max(0, lB.size() - lA.size()));
		return basesDifferent(getGraph().getK(), itA, itB);
	}
	/**
	 * Returns the number of bases difference between the two paths
	 * Bases are only compared up to the length of the shortest of the
	 * paths
	 * @param pathA
	 * @param pathB
	 * @return number of bases difference between the two paths
	 */
	protected int basesDifferent(Iterable<PN> pathA, Iterable<PN> pathB) {
		Iterator<Long> itA = PathNode.kmerIterator(pathA);
		Iterator<Long> itB = PathNode.kmerIterator(pathB);
		return basesDifferent(getGraph().getK(), itA, itB);
	}
	/**
	 * Returns the number of bases difference between the two paths
	 * Bases are only compared up to the length of the shortest of the
	 * paths
	 * @param k kmer length
	 * @param pathA
	 * @param pathB
	 * @return number of bases difference between the two paths
	 */
	private static int basesDifferent(int k, Iterator<Long> itA, Iterator<Long> itB) {
		int differences = 0;
		if (itA.hasNext() && itB.hasNext()) {
			differences = KmerEncodingHelper.basesDifference(k, itA.next(), itB.next());
			while (itA.hasNext() && itB.hasNext()) {
				if (!KmerEncodingHelper.lastBaseMatches(k, itA.next(), itB.next())) {
					differences++;
				}
			}
		}
		return differences;
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
					getGraph().getKmer(pathA.iterator().next().getFirst()).getSupportingEvidence().iterator().next().getBreakendSummary().toString(),
					new String(getGraph().getBaseCalls(Lists.newArrayList(PathNode.kmerIterator(pathA)))),
					new String(getGraph().getBaseCalls(Lists.newArrayList(PathNode.kmerIterator(pathB)))),
					new String(getGraph().getBaseCalls(Lists.newArrayList(PathNode.kmerIterator(ImmutableList.of(problematicNode)))))));
			}
			return false;
		}
		int weightA = PathNode.kmerTotalWeight(pathA);
		int weightB = PathNode.kmerTotalWeight(pathB);
		
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
	private PN getNodeAtDifferentPosition(Iterable<PN> pathA, Iterable<PN> pathB) {
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
			if (lengths.get(i) <= 0) throw new IllegalArgumentException("kmer path lengths must be greater than zero");
			sumlength += lengths.get(i);
		}
		if (sumlength != node.length()) throw new IllegalArgumentException("kmer path lengths must sum to node path length");
		// replace the source node with the split ones in the path graph
		List<PN> result = Lists.newArrayList();
		int offset = 0;
		for (int i = 0; i < lengths.size(); i++) {
			result.add(factory.createPathNode(node, offset, lengths.get(i), getGraph()));
			offset += lengths.get(i);
		}
		assert(offset == node.length());
		// now replace the node in the graph with the new nodes
		for (PN newNode : result) {
			addNode(newNode);
		}
		replaceIncomingEdges(node, result.get(0));
		for (int i = 0; i < result.size() - 1; i++) {
			nextPath(result.get(i)).add(result.get(i + 1));
			prevPath(result.get(i + 1)).add(result.get(i));
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
	 * A bubble is a de bruijn graph path diverges from a reference path
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
				return prevPath(node).size() == 1 && nextPath(node).size() == 1;
			}
		});
	}
	/**
	 * Traverses source graph kmers until a branch is found
	 * @param seed starting kmer
	 * @return unique path
	 */
	private LinkedList<Long> traverseBranchless(long seed) {
		LinkedList<Long> path = new LinkedList<Long>();
		Set<Long> visited = Sets.newHashSet();
		path.add(seed);
		visited.add(seed);
		for(List<Long> adj = getGraph().nextStates(path.getLast(), null, null); adj.size() == 1 && getGraph().prevStates(adj.get(0), null, null).size() <= 1; adj = getGraph().nextStates(path.getLast(), null, null)) {
			if (visited.contains(adj.get(0))) {
				// circular contig
				break;
			}
			path.addLast(adj.get(0));
			visited.add(adj.get(0));
		}
		for(List<Long> adj = getGraph().prevStates(path.getFirst(), null, null); adj.size() == 1 && getGraph().nextStates(adj.get(0), null, null).size() <= 1; adj = getGraph().prevStates(path.getFirst(), null, null)) {
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
		for (List<PN> nodeList = nextPath(path.getLast()); ; nodeList = nextPath(path.getLast())) {
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
		for (List<PN> nodeList = prevPath(path.getFirst()); ; nodeList = prevPath(path.getFirst())) {
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
	public PN getNodeContaining(long kmer) {
		for (PN node : paths) {
			if (node.contains(kmer)) return node;
		}
		return null;
	}
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	// Needs to be non-static to access graph.getK()
	public Ordering<PN> ByAverageKmerWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Doubles.compare(o1.getWeight() / (o1.length() + getGraph().getK() - 1), o2.getWeight() / (o2.length() + getGraph().getK() - 1));
		}
	}.reverse();
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	// Needs to be non-static as we can't have a generic type static property
	public Ordering<PN> ByPathTotalWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Ints.compare(o1.getWeight(), o2.getWeight());
		}
	}.reverse();
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	public Ordering<PN> ByMaxKmerWeightDesc = new Ordering<PN>() {
		public int compare(PN o1, PN o2) {
			return Ints.compare(o1.getMaxKmerWeight(), o2.getMaxKmerWeight());
		}
	}.reverse();
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("Path graph: %d node\n", paths.size()));
		for (PN x : paths) {
			sb.append(x.toString(getGraph()));
			sb.append(" <-{");
			for (PN s : prevPath(x)) sb.append(String.format("%d,", s.nodeId));
			sb.append("} ->{");
			for (PN s : nextPath(x)) sb.append(String.format("%d,", s.nodeId));
			sb.append("}\n");
		}
		return sb.toString();
	}
}
