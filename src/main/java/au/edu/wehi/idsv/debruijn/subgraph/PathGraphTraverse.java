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
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.graph.PathNode;
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
class PathGraphTraverse<T, PN extends PathNode<T>> {
	private static final Log log = Log.getInstance(PathGraphTraverse.class);
	private final Ordering<MemoizedNode<PN>> ByScore = new Ordering<MemoizedNode<PN>>() {
		@Override
		public int compare(MemoizedNode<PN> arg0, MemoizedNode<PN> arg1) {
			return Ints.compare(arg0.score, arg1.score);
		}
	};
	private final MemoizedNode<PN> NO_BEST_PATH = new MemoizedNode<PN>(null, null, Integer.MIN_VALUE, 0);
	private final DeBruijnPathGraph<T, PN> pg;
	private final int maxPathTraversalNodes;
	private final Function<PN, Integer> scoringFunction;
	private final Ordering<PN> order;
	private final boolean traverseForward;
	private final boolean referenceTraverse;
	private final int maxNextStates;
	private final int maxKmerPathLength;
	private int currentNextStates;
	private int nodeTraversals;
	private Map<PN, MemoizedNode<PN>> bestPath = Maps.newHashMap();
	private PriorityQueue<MemoizedNode<PN>> frontier = new PriorityQueue<MemoizedNode<PN>>(1024, ByScore);
	private MemoizedNode<PN> traversalBest;
	private Set<PN> terminalNodes;
	public int getNodesTraversed() { return nodeTraversals; }
	public PathGraphTraverse(
			final DeBruijnPathGraph<T, PN> pg,
			final int maxPathTraversalNodes,
			final Function<PN, Integer> scoringFunction,
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
	private static class MemoizedNode<PN> {
		public MemoizedNode(MemoizedNode<PN> prev, PN node, int score, int length) {
			this.prev = prev;
			this.node = node;
			this.score = score;
			this.length = length;
		}
		public final PN node;
		public final MemoizedNode<PN> prev;
		public final int score;
		public final int length;
	}
	public List<PN> traverse(final Iterable<PN> startingNodes) {
		return traverse(startingNodes, null);
	}
	public List<PN> traverse(final Iterable<PN> startingNodes, Collection<PN> endNodes) {
		if (endNodes == null) {
			terminalNodes = null;
		} else {
			if (endNodes instanceof Set<?>) {
				terminalNodes = (Set<PN>)endNodes;
			} else {
				terminalNodes = Sets.newHashSet(endNodes);
			}
		}
		currentNextStates = maxNextStates;
		nodeTraversals = 0;
		MemoizedNode<PN> globalBest = NO_BEST_PATH;
		for (PN starting : startingNodes) {
			MemoizedNode<PN> best;
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
		List<PN> bestPath = toList(globalBest, traverseForward); 
		return bestPath;
	}
	private MemoizedNode<PN> traverse(final PN startingNode) throws AlgorithmRuntimeSafetyLimitExceededException {
		init(startingNode);
		for (MemoizedNode<PN> node = popFronter(); node != null; node = popFronter()) {
			visitChildren(node);
		}
		return traversalBest;
	}
	private List<PN> toList(MemoizedNode<PN> node, boolean traverseForward) {
		if (node == NO_BEST_PATH) return null;
		List<PN> list = new ArrayList<PN>();
		while (node != null) {
			list.add(node.node);
			node = node.prev;
		}
		if (traverseForward) {
			list = Lists.reverse(list);
		}
		return list;
	}
	private void init(PN startingNode) {
		bestPath.clear();
		frontier.clear();
		traversalBest = NO_BEST_PATH;
		MemoizedNode<PN> start = new MemoizedNode<PN>(null, startingNode, SubgraphHelper.getScore(startingNode, scoringFunction, referenceTraverse, !referenceTraverse), startingNode.length());
		pushFrontier(start);
	}
	private void visitChildren(MemoizedNode<PN> node) throws AlgorithmRuntimeSafetyLimitExceededException {
		nodeTraversals++;
		if (nodeTraversals > maxPathTraversalNodes && currentNextStates > 1) {
			// maxNextStates = 1 is a greedy traversal - no need for a timeout as it is n nodes traversed even in a fully connected graph
			throw new AlgorithmRuntimeSafetyLimitExceededException();
		}
		assert(!Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS || sanityCheck());
		List<PN> nextStates = traverseForward ? pg.next(node.node) : pg.prev(node.node);
		Collections.sort(nextStates, order);
		int nextStateCount = 0;
		for (PN next : nextStates) {
			boolean isRef = pg.isReference(next);
			if ((referenceTraverse && !isRef) ||
				(!referenceTraverse && isRef)) {
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
	private boolean inPath(MemoizedNode<PN> path, PN node) {
		for (MemoizedNode<PN> mn = path; mn != null; mn = mn.prev) {
			if (node == mn.node) {
				return true;
			}
		}
		return false;
	}
	private boolean pushFrontier(MemoizedNode<PN> node) {
		MemoizedNode<PN> existing = bestPath.get(node.node);
		if (existing != null && existing.score >= node.score) {
			// our destination already has a better path to it - clearly not an optimal path 
			return false;
		}
		if (traversalBest.score < node.score) {
			if (terminalNodes == null || terminalNodes.contains(node.node)) {
				traversalBest = node;
			}
		}
		bestPath.put(node.node, node);
		frontier.add(node);
		return true;
	}
	private boolean pushFrontier(MemoizedNode<PN> node, PN next) {
		// TODO: calculate score only on bases within maxKmerPathLength
		int length = node.length + next.length();
		int score = node.score + scoringFunction.apply(next);
		return pushFrontier(new MemoizedNode<PN>(node, next, score, length));
	}
	private MemoizedNode<PN> popFronter() {
		while (!frontier.isEmpty()) {
			MemoizedNode<PN> node = frontier.poll();
			if (bestPath.get(node.node) == node && node.length < maxKmerPathLength) {
				return node;
			}
		}
		return null;
	}
	private boolean sanityCheck() {
		for (MemoizedNode<PN> n : frontier) {
			assert(terminalNodes != null || n.score <= traversalBest.score);
			assert(bestPath.containsKey(n.node));
			assert(bestPath.get(n.node).score >= n.score);
		}
		for (MemoizedNode<PN> n : bestPath.values()) {
			assert(bestPath.get(n.node).node == n.node);
		}
		return true;
	}
}

