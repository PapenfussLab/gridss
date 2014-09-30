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
	private int collapsePathTraversalCount;
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
	private boolean mergeWithNext(PN node) {
		if (nextPath(node).size() != 1) return false;
		PN next = nextPath(node).get(0); 
		if (next == node) return false; // don't collapse self
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
		sanityCheck();
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
	 */
	public void collapseSimilarPaths(int maxDifference, boolean bubblesOnly) {
		collapsePathTraversalCount = 0;
		try {
			// Keep collapsing until we can't anymore
			boolean collapsed = true;
			while (collapsed) {
				collapsed = false;
				for (PN start : paths) {
					// TODO: collapse the path with the weakest support first
					if (collapsePaths(maxDifference, bubblesOnly, start)) {
						collapsed = true;
						break;
					}
				}
			}
		} catch (AlgorithmRuntimeSafetyLimitExceededException e) {
			log.error(e);
		}
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
		if (++collapsePathTraversalCount >= COLLAPSE_PATH_MAX_TRAVERSAL) {
			throw new AlgorithmRuntimeSafetyLimitExceededException(String.format(
					"Path collapse not complete after %d node traversals. Aborting path collapse whilst processing \"%s\" and \"%s\".",
					collapsePathTraversalCount,
					new String(getGraph().getBaseCalls(Lists.newArrayList(PathNode.kmerIterator(pathA)))),
					new String(getGraph().getBaseCalls(Lists.newArrayList(PathNode.kmerIterator(pathB))))));
		}
		// paths have diverged too far
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
		assert(itA.hasNext());
		assert(itB.hasNext());
		int differences = KmerEncodingHelper.basesDifference(getGraph().getK(), itA.next(), itB.next());
		while (itA.hasNext() && itB.hasNext()) {
			if (!KmerEncodingHelper.lastBaseMatches(getGraph().getK(), itA.next(), itB.next())) {
				differences++;
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
		int weightA = getWeight(pathA);
		int weightB = getWeight(pathB);
		
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
				// contains both: just remove
				list.remove(toReplace);
			} else {
				// only contains toReplace, just overwrite with the new value
				list.set(list.indexOf(toReplace), replaceWith);
			}
		}
	}
	/**
	 * Gets the total weight of all kmers on the given path
	 * @param path path
	 * @return total weight
	 */
	public int getWeight(Iterable<PN> path) {
		int weight = 0;
		for (PN n : path) {
			weight += n.getWeight();
		}
		return weight;
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
