package au.edu.wehi.idsv.debruijn;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;

import au.edu.wehi.idsv.graph.WeightedSequenceGraphNodeUtil;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.NodeTraversalTracker;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class HighestWeightSimilarPath<PN extends DeBruijnSequenceGraphNode> {
	private List<PN> bestPath = null;
	private int bestWeight = Integer.MIN_VALUE;
	private int maxDifference;
	private PN leaf;
	private PN anchor;
	private boolean traverseForward;
	private DeBruijnGraph<PN> graph;
	private NodeTraversalTracker<PN> tracker;
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
	public HighestWeightSimilarPath(
			int maxDifference,
			PN leaf,
			boolean traverseForward,
			DeBruijnGraph<PN> graph,
			NodeTraversalTracker<PN> tracker) {
		this.bestPath = null;
		this.maxDifference = maxDifference;
		this.leaf = leaf;
		List<PN> anchorList = traverseForward ? graph.prev(leaf) : graph.next(leaf);
		assert(anchorList.size() == 1); // must be leaf node
		this.anchor = anchorList.get(0);
		this.traverseForward = traverseForward;
		this.graph = graph;
		this.tracker = tracker;
	}
	public List<PN> find() {
		recursiveFind(new ArrayDeque<PN>(), 0, 0);
		return bestPath;
	}
	private void recursiveFind(
			Deque<PN> currentPath,
			int currentLength,
			int currentDifferences) {
		if (currentDifferences > maxDifference) {
			return;
		}
		List<PN> nextList = traverseForward ? graph.next((currentPath.isEmpty() ? anchor : currentPath.getLast())) : graph.prev((currentPath.isEmpty() ? anchor : currentPath.getFirst()));
		if (currentLength >= leaf.length()) {
			assert(currentLength == WeightedSequenceGraphNodeUtil.nodeLength(currentPath));
			int leafLength = leaf.length();
			int currentSharedLength = Math.min(currentLength, leafLength);
			int pathWeight = WeightedSequenceGraphNodeUtil.totalWeight(currentPath, traverseForward ? 0 : currentLength - currentSharedLength, currentSharedLength);
			if (pathWeight > bestWeight) {
				// this path has a higher total weight over the kmers shared with the leaf 
				bestPath = Lists.newArrayList(currentPath);
				bestWeight = pathWeight;
			} else {
				return;
			}
		}
		if (nextList.size() == 0) {
			// Don't merge a leaf into shorter paths
			// Doing so messes up consensus kmer sequence - the primary sequence
			// would no longer be a valid kmer path. We also can't try to fix the
			// consensus sequence as the kmer we replace the inconsistent kmers
			// with could already be part of our graph elsewhere
			// or, even worse, part of a different subgraph in the kmer graph! 
			return;
		}
		for (PN next : nextList) {
			if (next != leaf && !currentPath.contains(next)) {
				if (traverseForward) {
					currentPath.addLast(next);
				} else {
					currentPath.addFirst(next);
				}
				if (tracker != null) {
					tracker.track(next);
				}
				recursiveFind(
						currentPath,
						currentLength + next.length(),
						traverseForward ? DeBruijnSequenceGraphNodeUtil.basesDifferent(graph.getK(), ImmutableList.of(leaf), currentPath) : DeBruijnSequenceGraphNodeUtil.reverseBasesDifferent(graph.getK(), ImmutableList.of(leaf), currentPath)
						);
				if (traverseForward) {
					currentPath.removeLast();
				} else {
					currentPath.removeFirst();
				}
			}
		}
	}
}
