package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;

import com.google.common.base.Function;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;
/**
 * Dynamic programming traversal of subgraph of path graph
 * 
 * @author cameron.d
 *
 */
class PathGraphTraverse {
	private static final Log log = Log.getInstance(PathGraphTraverse.class);
	private static final MemoizedNode noBestPath = new MemoizedNode(null, null, Integer.MIN_VALUE, 0);
	private final PathGraph pg;
	private final int maxPathTraversalNodes;
	private final Function<SubgraphPathNode, Integer> scoringFunction;
	private final Ordering<SubgraphPathNode> order;
	private final boolean traverseForward;
	private final boolean referenceTraverse;
	private final int maxNextStates;
	private final int maxKmerPathLength;
	private int currentNextStates;
	private int nodeTraversals;
	private Map<SubgraphPathNode, MemoizedNode> bestPath = Maps.newHashMap();
	private PriorityQueue<MemoizedNode> frontier = new PriorityQueue<MemoizedNode>(1024, ByScore);
	private MemoizedNode traveralBest;
	private Set<SubgraphPathNode> terminalNodes;
	public int getNodesTraversed() { return nodeTraversals; }
	public PathGraphTraverse(
			final PathGraph pg,
			final int maxPathTraversalNodes,
			final Function<SubgraphPathNode, Integer> scoringFunction,
			final boolean traverseForward,
			final boolean referenceTraverse,
			final int maxNextStates,
			final int maxKmerPathLength) {
		this.pg = pg;
		this.maxPathTraversalNodes = maxPathTraversalNodes;
		this.scoringFunction = scoringFunction;
		this.order = SubgraphHelper.getScoreOrderDesc(scoringFunction);
		this.traverseForward = traverseForward;
		this.referenceTraverse = referenceTraverse;
		this.maxNextStates = maxNextStates;
		this.maxKmerPathLength = maxKmerPathLength;
	}
	private static final Ordering<MemoizedNode> ByScore = new Ordering<MemoizedNode>() {
		@Override
		public int compare(MemoizedNode arg0, MemoizedNode arg1) {
			return Ints.compare(arg0.score, arg1.score);
		}
	};
	private static final Ordering<MemoizedNode> ByLength = new Ordering<MemoizedNode>() {
		@Override
		public int compare(MemoizedNode arg0, MemoizedNode arg1) {
			return Ints.compare(arg0.length, arg1.length);
		}
	};
	private static class MemoizedNode {
		public MemoizedNode(MemoizedNode prev, SubgraphPathNode node, int score, int length) {
			this.prev = prev;
			this.node = node;
			this.score = score;
			this.length = length;
		}
		public final SubgraphPathNode node;
		public final MemoizedNode prev;
		public final int score;
		public final int length;
	}
	public List<SubgraphPathNode> traverse(final Iterable<SubgraphPathNode> startingNodes) {
		return traverse(startingNodes, null);
	}
	public List<SubgraphPathNode> traverse(final Iterable<SubgraphPathNode> startingNodes, Collection<SubgraphPathNode> endNodes) {
		if (endNodes == null) {
			terminalNodes = null;
		} else {
			if (endNodes instanceof Set<?>) {
				terminalNodes = (Set<SubgraphPathNode>)endNodes;
			} else {
				terminalNodes = Sets.newHashSet(endNodes);
			}
		}
		currentNextStates = maxNextStates;
		nodeTraversals = 0;
		MemoizedNode globalBest = noBestPath;
		for (SubgraphPathNode starting : startingNodes) {
			MemoizedNode best;
			try {
				best = traverse(starting);
			} catch (AlgorithmRuntimeSafetyLimitExceededException e) {
				currentNextStates = 1;
				log.debug(String.format("Reached path traversal timeout. Switching to greedy traversal."));
				try {
					best = traverse(starting);
				} catch (AlgorithmRuntimeSafetyLimitExceededException e1) {
					log.error("Sanity check failure: timeout hit again when performing greedy traversal");
					throw new RuntimeException(e);
				}
			}
			if (best.score > globalBest.score) {
				globalBest = best;
			}
		}
		List<SubgraphPathNode> bestPath = toList(globalBest, traverseForward); 
		return bestPath;
	}
	private MemoizedNode traverse(final SubgraphPathNode startingNode) throws AlgorithmRuntimeSafetyLimitExceededException {
		init(startingNode);
		for (MemoizedNode node = popFronter(); node != null; node = popFronter()) {
			visit(node);
		}
		return traveralBest;
	}
	private static List<SubgraphPathNode> toList(MemoizedNode node, boolean traverseForward) {
		if (node == noBestPath) return null;
		List<SubgraphPathNode> list = new ArrayList<SubgraphPathNode>();
		while (node != null) {
			list.add(node.node);
			node = node.prev;
		}
		if (traverseForward) {
			list = Lists.reverse(list);
		}
		return list;
	}
	private void init(SubgraphPathNode startingNode) {
		bestPath.clear();
		frontier.clear();
		traveralBest = noBestPath;
		MemoizedNode start = new MemoizedNode(null, startingNode, SubgraphHelper.getScore(startingNode, scoringFunction, referenceTraverse, !referenceTraverse), startingNode.length());
		bestPath.put(startingNode, start);
		frontier.add(start);
	}
	private void visit(MemoizedNode node) throws AlgorithmRuntimeSafetyLimitExceededException {
		nodeTraversals++;
		if (nodeTraversals > maxPathTraversalNodes && currentNextStates > 1) {
			// maxNextStates = 1 is a greedy traversal - no need for a timeout as it is n nodes traversed even in a fully connected graph
			throw new AlgorithmRuntimeSafetyLimitExceededException();
		}
		assert(!Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS || sanityCheck());
		if (traveralBest.score < node.score) {
			if (terminalNodes == null || terminalNodes.contains(node.node)) {
				traveralBest = node;
			}
		}
		List<SubgraphPathNode> nextStates = traverseForward ? pg.nextPath(node.node) : pg.prevPath(node.node);
		Collections.sort(nextStates, order);
		int nextStateCount = 0;
		for (SubgraphPathNode next : nextStates) {
			if ((referenceTraverse && !next.containsReferenceKmer()) ||
				(!referenceTraverse && !next.containsNonReferenceKmer())) {
				// don't want to transition between reference and non-reference
				continue;
			}
			// path must be simple path
			if (inPath(node, next)) continue;
			
			// ok, we have somewhere to go!
			if (pushFrontier(node, next)) {
				nextStateCount++;
				if (nextStateCount >= currentNextStates) {
					// out we go: we've tried enough paths
					break;
				}
			}
		}
	}
	/**
	 * Determines whether the given node is in the given path
	 * @param path
	 * @param node
	 * @return
	 */
	private boolean inPath(MemoizedNode path, SubgraphPathNode node) {
		for (MemoizedNode mn = path; mn != null; mn = mn.prev) {
			if (node == mn.node) {
				return true;
			}
		}
		return false;
	}
	private boolean pushFrontier(MemoizedNode node, SubgraphPathNode next) {
		// TODO: calculate score only on bases within maxKmerPathLength
		int length = node.length + next.length();
		int score = node.score + scoringFunction.apply(next);
		MemoizedNode existing = bestPath.get(next);
		if (existing != null && existing.score >= score) {
			// our destination already has a better path to it - clearly not an optimal path 
			return false;
		}
		MemoizedNode newNode = new MemoizedNode(node, next, score, length);
		bestPath.put(next, newNode);
		frontier.add(newNode);
		return true;
	}
	private MemoizedNode popFronter() {
		while (!frontier.isEmpty()) {
			MemoizedNode node = frontier.poll();
			if (bestPath.get(node.node) == node && node.length < maxKmerPathLength) {
				return node;
			}
		}
		return null;
	}
	private boolean sanityCheck() {
		for (MemoizedNode n : frontier) {
			assert(terminalNodes != null || n.score <= traveralBest.score);
			assert(bestPath.containsKey(n.node));
			assert(bestPath.get(n.node).score >= n.score);
		}
		for (MemoizedNode n : bestPath.values()) {
			assert(bestPath.get(n.node).node == n.node);
		}
		return true;
	}
}

