package au.edu.wehi.idsv.debruijn;

import java.util.Collection;
import java.util.Deque;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.graph.PathGraph;
import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.graph.PathNodeFactory;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;

import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Queues;

public class DeBruijnPathGraph<T, PN extends PathNode<T>> extends PathGraph<T, PN> implements DeBruijnGraph<PN> {
	private int nodeTaversalCount;
	public DeBruijnPathGraph(DeBruijnGraph<T> graph, Collection<T> seeds, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		super(graph, seeds, factory, tracker);
	}
	public DeBruijnPathGraph(DeBruijnGraph<T> graph, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		super(graph, factory, tracker);
	}
	@Override
	public DeBruijnGraph<T> getGraph() {
		return (DeBruijnGraph<T>)super.getGraph();
	}
	@Override
	public boolean isReference(PN node) {
		// TODO: cache this value in the node itself
		return isReference(node.allNodes());
	}
	private boolean isReference(Iterable<T> nodeSet) {
		DeBruijnGraph<T> g = getGraph();
		for (T n : nodeSet) {
			if (g.isReference(n)) {
				return true;
			}
		}
		return false;
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
	public int collapseSimilarPaths(int maxDifference, boolean bubblesOnly) throws AlgorithmRuntimeSafetyLimitExceededException {
		nodeTaversalCount = 0;
		int collapseCount = 0;
		// Keep collapsing until we can't anymore
		boolean collapsed = true;
		while (collapsed) {
			collapsed = false;
			for (PN start : getPaths()) {
				// TODO: collapse the path with the weakest support first
				if (collapsePaths(maxDifference, bubblesOnly, start)) {
					collapseCount++;
					collapsed = true;
					break;
				}
			}
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
	public int collapseLeaves(int maxDifference) throws AlgorithmRuntimeSafetyLimitExceededException {
		nodeTaversalCount = 0;
		int collapseCount = 0;
		int collapsedThisRound;
		do {
			collapsedThisRound = 0;
			PriorityQueue<PN> ordering = new PriorityQueue<PN>(getPathCount(), ByPathTotalWeight);
			for (PN p : getPaths()) ordering.add(p); //ordering.addAll(getPaths()); // no PriorityQueue.addAll(Iterable) method :(
			// unlike collapsing paths, collapsing leaves is a local
			// change and does not invalidate any other paths
			for (PN leaf : ordering) {
				List<PN> nexts = next(leaf);
				List<PN> prevs = prev(leaf);
				if (nexts.size() == 0 && prevs.size() == 1) {
					if (collapseLeaf(maxDifference, leaf, prevs.get(0), true)) {
						collapseCount++;
						collapsedThisRound++;
					}
				} else if (nexts.size() == 1 && prevs.size() == 0) {
					if (collapseLeaf(maxDifference, leaf, nexts.get(0), false)) {
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
		return collapseCount;
	}
	/**
	 * Collapses similar paths and leaf branches
	 * @param maxBaseMismatchForCollapse
	 * @param bubblesOnly collapse bubbles only
	 * @return total collapses performed
	 */
	public int collapse(int maxBaseMismatchForCollapse, boolean bubblesOnly) throws AlgorithmRuntimeSafetyLimitExceededException {
		int nodesBefore = getPathCount();
		int totalPathsCollapsed = 0;
		int totalLeavesCollapsed = 0;
		int collapseIterations = 0;
		totalPathsCollapsed = collapseSimilarPaths(maxBaseMismatchForCollapse, true);
		int collapseLeafCount = 0, collapsePathCount = 0;
		do {
			collapseIterations++;
			if (!bubblesOnly) {
				collapseLeafCount = collapseLeaves(maxBaseMismatchForCollapse);
				totalLeavesCollapsed += collapseLeafCount;
			}
			collapsePathCount = collapseSimilarPaths(maxBaseMismatchForCollapse, bubblesOnly);
			totalPathsCollapsed += collapsePathCount;
		} while (collapseLeafCount + collapsePathCount > 0);
		tracker.collapse(collapseIterations, totalPathsCollapsed, totalLeavesCollapsed, nodeTaversalCount, nodesBefore - getPathCount());
		return totalPathsCollapsed + totalLeavesCollapsed;
	}
	/**
	 * 
	 * @param maxDifference maximum bases different
	 * @param leaf leaf node
	 * @param anchor anchor to start collapse
	 * @param forward anchor to leaf is in the direction of kmer traversal 
	 * @return true if the leaf was collapsed, false if no similar path found
	 * @throws AlgorithmRuntimeSafetyLimitExceededException 
	 */
	private boolean collapseLeaf(int maxDifference, PN leaf, PN anchor, boolean forward) throws AlgorithmRuntimeSafetyLimitExceededException {
		// find path from anchor that is within maxDifference from leaf
		List<PN> path = findHighestWeightSimilarPath(ImmutableList.<PN>of(), Queues.<PN>newArrayDeque(), 0, 0, maxDifference, leaf, anchor, forward);
		if (path == null || path.isEmpty()) return false;
		// compare weights of shared kmers
		int pathLength = PathNode.nodeLength(path);
		int sharedLength = Math.min(leaf.length(), pathLength);
		if (PathNode.totalWeight(ImmutableList.of(leaf), forward ? 0 : leaf.length() - sharedLength, sharedLength, getGraph()) >=
			PathNode.totalWeight(path, forward ? 0 : pathLength - sharedLength, sharedLength, getGraph())) {
			return false;
		}
		collapseLeafInto(leaf, path, forward);
		assert(sanityCheck());
		return true;
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
		if (++nodeTaversalCount >= Defaults.COLLAPSE_PATH_MAX_TRAVERSAL) {
			throw new AlgorithmRuntimeSafetyLimitExceededException(String.format(
					"Leaf collapse not complete after %d node traversals. Aborting path collapse whilst processing \"%s\".",
					nodeTaversalCount,
					toString(ImmutableList.of(leaf))));
		}
		if (currentDifferences > maxDifference) {
			return bestPath;
		}
		List<PN> nextList = traverseForward ? next((currentPath.isEmpty() ? anchor : currentPath.getLast())) : prev((currentPath.isEmpty() ? anchor : currentPath.getFirst()));
		if (currentLength >= leaf.length()) {
			assert(currentLength == PathNode.nodeLength(currentPath));
			int leafLength = leaf.length();
			int bestLength = PathNode.nodeLength(bestPath);
			int currentSharedLength = Math.min(currentLength, leafLength);
			int bestSharedLength = Math.min(bestLength, leafLength);
			if (PathNode.totalWeight(currentPath, traverseForward ? 0 : currentLength - currentSharedLength, currentSharedLength, getGraph()) >
				PathNode.totalWeight(bestPath, traverseForward ? 0 : bestLength - bestSharedLength, bestSharedLength, getGraph())) {
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
		List<PN> next = next(start);
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
		if (++nodeTaversalCount >= Defaults.COLLAPSE_PATH_MAX_TRAVERSAL) {
			throw new AlgorithmRuntimeSafetyLimitExceededException(String.format(
					"Path collapse not complete after %d node traversals. Aborting path collapse whilst processing \"%s\" and \"%s\".",
					nodeTaversalCount,
					toString(pathA),
					toString(pathB)));
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
			for (PN nextA : next(pathA.getLast())) {
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
			for (PN nextB : next(pathB.getLast())) {
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
	 * Returns the number of bases difference between the two paths when
	 * traversing the paths backward.
	 * Bases are only compared up to the length of the shortest of the
	 * paths
	 * @param pathA
	 * @param pathB
	 * @return number of bases difference between the two paths
	 */
	public int reverseBasesDifferent(Iterable<PN> pathA, Iterable<PN> pathB) {
		int lengthA = PathNode.nodeLength(pathA);
		int lengthB = PathNode.nodeLength(pathB);
		int skipCountA = Math.max(0, lengthA - lengthB);
		int skipCountB = Math.max(0, lengthB - lengthA);
		// skip initial bases of the longer path
		return KmerEncodingHelper.totalBaseDifference(getGraph(),
				PathNode.nodeIterator(pathA.iterator(), skipCountA),
				PathNode.nodeIterator(pathB.iterator(), skipCountB));
	}
	/**
	 * Returns the number of bases difference between the two paths
	 * Bases are only compared up to the length of the shortest of the
	 * paths
	 * @param pathA
	 * @param pathB
	 * @return number of bases difference between the two paths
	 */
	public int basesDifferent(Iterable<PN> pathA, Iterable<PN> pathB) {
		return KmerEncodingHelper.totalBaseDifference(getGraph(),
				PathNode.nodeIterator(pathA.iterator()),
				PathNode.nodeIterator(pathB.iterator()));
	}
	private int referencePathsSplit = 0;
	public int getReferencePathsSplit() { return referencePathsSplit; }
	/**
	 * Breaks paths into paths in which all kmers are either reference or non-reference kmers 
	 */
	public void splitOutReferencePaths() {
		// get full list of splits required upfront (since we'll be changing the node set as we iterate) 
		List<PN> toSplit = Lists.newArrayList(Iterables.filter(getPaths(), new Predicate<PN>() {
			public boolean apply(PN arg) {
				boolean containsReference = false;
				boolean containsNonReference = false;
				for (List<T> nodeSet : arg.getPathAllNodes()) {
					boolean isRef = isReference(nodeSet);
					containsReference |= isRef;
					containsNonReference |= !isRef;
					if (containsReference && containsNonReference) {
						return true;
					}
				}
				return false;
			}
		}));
		for (PN n : toSplit) {
			List<Integer> lengths = Lists.newArrayList();
			List<List<T>> path = n.getPathAllNodes();
			List<Boolean> ref = Lists.newArrayList();
			int currentLength = 0;
			boolean currentIsReference = isReference(path.get(0)); 
			for (List<T> kmers : path) {
				boolean isRef = isReference(kmers);
				if (currentIsReference == isRef) {
					currentLength++;
				} else {
					currentIsReference = isRef;
					lengths.add(currentLength);
					currentLength = 1;
					ref.add(isRef);
				}
			}
			ref.add(currentIsReference);
			lengths.add(currentLength);
			split(n, lengths);
			referencePathsSplit += lengths.size() - 1;
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(assertReferenceKmersSplit());
		}
	}
	private boolean assertReferenceKmersSplit() {
		for (PN pn : getPaths()) {
			boolean shouldBeReference = isReference(pn);
			for (List<T> list : pn.getPathAllNodes()) {
				int referenceNodeCount = 0;
				for (T n : list) {
					if (getGraph().isReference(n)) {
						referenceNodeCount++;
					}
				}
				boolean containsReference = referenceNodeCount > 0;
				assert(containsReference == shouldBeReference);
			}
		}
		return true;
	}
	@Override
	public long getKmer(PN node) {
		throw new UnsupportedOperationException("path nodes contain multiple kmers");
	}
	@Override
	public int getK() {
		return getGraph().getK();
	}
}
