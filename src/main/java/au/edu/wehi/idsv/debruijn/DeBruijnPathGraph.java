package au.edu.wehi.idsv.debruijn;

import java.util.Collection;
import java.util.Deque;
import java.util.List;
import java.util.PriorityQueue;

import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Queues;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.graph.PathGraph;
import au.edu.wehi.idsv.graph.PathNodeFactory;
import au.edu.wehi.idsv.graph.WeightedSequenceGraphNodeUtil;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;
import au.edu.wehi.idsv.visualisation.TimeoutNodeTraversalTracker;

public class DeBruijnPathGraph<T, PN extends DeBruijnPathNode<T>> extends PathGraph<T, PN> implements DeBruijnGraph<PN> {
	private TimeoutNodeTraversalTracker<PN> collapseTimeout;
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
		return node.isReference();
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
	public int collapseSimilarPaths(int maxDifference, boolean bubblesOnly, int maxTraversalNodes) throws AlgorithmRuntimeSafetyLimitExceededException {
		collapseTimeout = new TimeoutNodeTraversalTracker<PN>(maxTraversalNodes);
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
	public int collapseLeaves(int maxDifference, int maxTraversalNodes) throws AlgorithmRuntimeSafetyLimitExceededException {
		collapseTimeout = new TimeoutNodeTraversalTracker<PN>(maxTraversalNodes);
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
	public int collapse(int maxBaseMismatchForCollapse, boolean bubblesOnly, int maxTraversalNodes) throws AlgorithmRuntimeSafetyLimitExceededException {
		int nodesBefore = getPathCount();
		int totalPathsCollapsed = 0;
		int totalLeavesCollapsed = 0;
		int collapseIterations = 0;
		totalPathsCollapsed = collapseSimilarPaths(maxBaseMismatchForCollapse, true, maxTraversalNodes);
		int collapseLeafCount = 0, collapsePathCount = 0;
		do {
			collapseIterations++;
			if (!bubblesOnly) {
				collapseLeafCount = collapseLeaves(maxBaseMismatchForCollapse, maxTraversalNodes);
				totalLeavesCollapsed += collapseLeafCount;
			}
			collapsePathCount = collapseSimilarPaths(maxBaseMismatchForCollapse, bubblesOnly, maxTraversalNodes);
			totalPathsCollapsed += collapsePathCount;
		} while (collapseLeafCount + collapsePathCount > 0);
		tracker.collapse(collapseIterations, totalPathsCollapsed, totalLeavesCollapsed, collapseTimeout.count(), nodesBefore - getPathCount());
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
		List<PN> path = new HighestWeightSimilarPath<PN>(maxDifference, leaf, forward, this, collapseTimeout).find();
		if (path == null || path.isEmpty()) return false;
		// compare weights of shared kmers
		int pathLength = WeightedSequenceGraphNodeUtil.nodeLength(path);
		int sharedLength = Math.min(leaf.length(), pathLength);
		if (WeightedSequenceGraphNodeUtil.totalWeight(ImmutableList.of(leaf), forward ? 0 : leaf.length() - sharedLength, sharedLength) >=
				WeightedSequenceGraphNodeUtil.totalWeight(path, forward ? 0 : pathLength - sharedLength, sharedLength)) {
			return false;
		}
		collapseLeafInto(leaf, path, forward);
		assert(sanityCheck());
		return true;
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
		// paths have diverged too far
		// Optimisation: cache bases compared and difference so far so only 2 PNs need to be compared
		if (DeBruijnSequenceGraphNodeUtil.basesDifferent(getK(), pathA, pathB) > differencesAllowed) return false;
		
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
					collapseTimeout.track(nextA);
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
					collapseTimeout.track(nextB);
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
		if (Defaults.SANITY_CHECK_DE_BRUIJN) {
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
