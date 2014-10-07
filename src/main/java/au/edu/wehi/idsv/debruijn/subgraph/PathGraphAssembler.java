package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.util.Log;

import java.util.ArrayDeque;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.Queue;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.visualisation.DeBruijnPathGraphExporter;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

/**
 * Reduces the de bruijn graph to a set of paths before assembling
 * @author Daniel Cameron
 *
 */
public class PathGraphAssembler extends PathGraph {
	private final AssemblyParameters parameters;
	private final List<Set<SubgraphPathNode>> subgraphs = Lists.newArrayList();
	private final List<Set<SubgraphPathNode>> startingPaths = Lists.newArrayList();
	public PathGraphAssembler(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, AssemblyParameters parameters, long seed) {
		super(graph, seed);
		this.parameters = parameters;
	}
	public List<LinkedList<Long>> assembleContigs() {
		return assembleContigs(null);
	}
	public List<LinkedList<Long>> assembleContigs(DeBruijnPathGraphExporter<DeBruijnSubgraphNode, SubgraphPathNode> graphExporter) {
		splitOutReferencePaths();
		calcNonReferenceSubgraphs();
		return assemblyNonReferenceContigs(graphExporter);
	}
	private List<LinkedList<Long>> assemblyNonReferenceContigs(DeBruijnPathGraphExporter<DeBruijnSubgraphNode, SubgraphPathNode> graphExporter) {
		final Function<SubgraphPathNode, Integer> scoringFunction = SCORE_TOTAL_KMER;
		List<List<SubgraphPathNode>> result = Lists.newArrayList();
		for (int i = 0; i < subgraphs.size(); i++) {
			List<SubgraphPathNode> contig = assembleSubgraph(subgraphs.get(i), startingPaths.get(i), scoringFunction);
			result.add(contig);
		}
		// return the best assembly first
		Collections.sort(result, new Ordering<List<SubgraphPathNode>>() {
			@Override
			public int compare(List<SubgraphPathNode> arg0, List<SubgraphPathNode> arg1) {
				return Ints.compare(
						getScore(arg1, scoringFunction, false, true),
						getScore(arg0, scoringFunction, false, true));
			}});
		if (graphExporter != null) {
			graphExporter.snapshot(this);
			graphExporter.annotateSubgraphs(subgraphs);
			graphExporter.annotateStartingPaths(startingPaths);
			graphExporter.contigs(result);
		}
		return Lists.newArrayList(Iterables.transform(result, new Function<List<SubgraphPathNode>, LinkedList<Long>>() {
			@Override
			public LinkedList<Long> apply(List<SubgraphPathNode> arg0) {
				return Lists.newLinkedList(Lists.newArrayList(SubgraphPathNode.kmerIterator(arg0)));
			}
		}));
	}
	/**
	 * Separates non-reference kmers of the subgraph into disjoint non-reference subgraphs 
	 */
	private void calcNonReferenceSubgraphs() {
		Set<SubgraphPathNode> toProcess = Sets.newHashSet(getPaths());
		while (!toProcess.isEmpty()) {
			SubgraphPathNode node = toProcess.iterator().next();
			toProcess.remove(node);
			if (node.containsNonReferenceKmer()) {
				newNonreferenceSubgraph(node);
				toProcess.removeAll(subgraphs.get(subgraphs.size() - 1));
			}
		}
		if (DeBruijnGraphBase.PERFORM_EXPENSIVE_SANITY_CHECKS) {
			assert(sanityCheckSubgraphs());
		}
	}
	private void newNonreferenceSubgraph(SubgraphPathNode seed) {
		Set<SubgraphPathNode> nodes = Sets.newHashSet();
		Set<SubgraphPathNode> startingNodes = Sets.newHashSet();
		Queue<SubgraphPathNode> frontier = new ArrayDeque<SubgraphPathNode>();
		nodes.add(seed);
		frontier.add(seed);
		while (!frontier.isEmpty()) {
			SubgraphPathNode node = frontier.poll();
			for (SubgraphPathNode prev : prevPath(node)) {
				if (prev.containsReferenceKmer()) {
					startingNodes.add(node);
				}
				if (prev.containsNonReferenceKmer() && !nodes.contains(prev)) {
					nodes.add(prev);
					frontier.add(prev);
				}
			}
			for (SubgraphPathNode next : nextPath(node)) {
				if (next.containsNonReferenceKmer() && !nodes.contains(next)) {
					nodes.add(next);
					frontier.add(next);
				}
			}
		}
		if (startingNodes.isEmpty()) { 
			// no reference anchors at all - try starting anywhere
			// TODO: reduce # starting positions by restricting to either
			// 1) no incoming edges (in the nodes set)
			// 2) nodes forming part of a cycle 
			startingNodes = nodes;
		}
		subgraphs.add(nodes);
		startingPaths.add(startingNodes);
	}
	private List<SubgraphPathNode> assembleSubgraph(
			Set<SubgraphPathNode> nodes,
			Set<SubgraphPathNode> startingNodes,
			Function<SubgraphPathNode, Integer> scoringFunction) {
		int branchingFactor = parameters.subgraphAssemblyTraversalMaximumBranchingFactor;
		List<SubgraphPathNode> best = null;
		int bestScore = Integer.MIN_VALUE;
		for (SubgraphPathNode startingNode : startingNodes) {
			List<SubgraphPathNode> current = memoizedTraverse(startingNode, scoringFunction, true, false, branchingFactor);
			int currentScore = getScore(current, scoringFunction, false, true);
			if (currentScore >= bestScore) {
				bestScore = currentScore;
				best = current;
			}
		}
		List<SubgraphPathNode> anchor = memoizedTraverse(best.get(0), scoringFunction, false, true, branchingFactor);
		anchor.remove(anchor.size() - 1); // don't repeat the non-reference node we used as the anchor
		return Lists.newArrayList(Iterables.concat(anchor, best));
	}
	private int getScore(
			List<SubgraphPathNode> current,
			Function<SubgraphPathNode,
			Integer> perNodeScoringFunction,
			boolean scoreReference,
			boolean scoreNonReference) {
		int sum = 0;
		for (SubgraphPathNode n : current) {
			if ((n.containsNonReferenceKmer() && scoreNonReference) ||
				 (n.containsReferenceKmer() && scoreReference)) {
				sum += perNodeScoringFunction.apply(n);
			}
		}
		return sum;
	}
	private static Ordering<SubgraphPathNode> getScoreOrderDesc(final Function<SubgraphPathNode, Integer> scoringFunction) {
		return new Ordering<SubgraphPathNode>() {
			@Override
			public int compare(SubgraphPathNode arg0, SubgraphPathNode arg1) {
				return Ints.compare(scoringFunction.apply(arg1), scoringFunction.apply(arg0));
			}};
	}
	private List<SubgraphPathNode> memoizedTraverse(
			final SubgraphPathNode startingNode,
			final Function<SubgraphPathNode, Integer> scoringFunction,
			final boolean traverseForward,
			final boolean referenceTraverse,
			final int maxNextStates) {
		Ordering<SubgraphPathNode> order = getScoreOrderDesc(scoringFunction);
			
		Map<SubgraphPathNode, Integer> bestScore = Maps.newHashMap();
		Map<SubgraphPathNode, List<SubgraphPathNode>> bestPath = Maps.newHashMap();
		NavigableMap<Integer, Set<SubgraphPathNode>> frontier = Maps.newTreeMap();
		int bestFinalScore = Integer.MIN_VALUE;
		List<SubgraphPathNode> bestFinalPath = null;
		List<SubgraphPathNode> startingPath = Lists.newArrayList(startingNode); 
		int startingScore = getScore(startingPath, scoringFunction, referenceTraverse, !referenceTraverse);
		bestScore.put(startingNode, startingScore);
		bestPath.put(startingNode, startingPath);
		pushFronter(frontier, startingNode, startingScore, null);
		
		while (!frontier.isEmpty()) {
			if (DeBruijnGraphBase.PERFORM_EXPENSIVE_SANITY_CHECKS) {
				assert(sanityCheckMemoization(bestScore, bestPath, frontier, scoringFunction, referenceTraverse, !referenceTraverse));
			}
			SubgraphPathNode node = popFronter(frontier);
			int score = bestScore.get(node);
			List<SubgraphPathNode> nodePath = bestPath.get(node);
			List<SubgraphPathNode> nextStates = traverseForward ? nextPath(node) : prevPath(node);
			Collections.sort(nextStates, order);
			int nextStateCount = 0;
			for (SubgraphPathNode next : nextStates) {
				if ((referenceTraverse && !next.containsReferenceKmer()) ||
					(!referenceTraverse && !next.containsNonReferenceKmer())) {
					// don't want to transition between reference and non-reference
					continue;
				}
				// path must be simple path 
				if (nodePath.contains(next)) continue;
				
				// ok, we have somewhere to go!
				int nextScore = score + scoringFunction.apply(next);
				if (bestScore.containsKey(next) && (int)bestScore.get(next) >= nextScore) {
					// our destination already has a better path to it - clearly not an optimal path
					continue;
				}
				nextStateCount++;
				if (nextStateCount > maxNextStates) {
					// out we go: we've tried enough times
					continue;
				}
				// memoize our path to our next node
				List<SubgraphPathNode> nextPath = Lists.newArrayListWithCapacity(nodePath.size() + 1);
				nextPath.addAll(nodePath);
				nextPath.add(next);
				pushFronter(frontier, next, nextScore, bestScore.containsKey(next) ? bestScore.get(next) : null);
				bestScore.put(next, nextScore);
				bestPath.put(next, nextPath);
			}
			if (nextStateCount == 0 && score > bestFinalScore) {
				bestFinalScore = score;
				bestFinalPath = nodePath;
			}
		}
		if (!traverseForward) {
			bestFinalPath = Lists.reverse(bestFinalPath);
		}
		return bestFinalPath;
	}
	private static void pushFronter(NavigableMap<Integer, Set<SubgraphPathNode>> frontier, SubgraphPathNode node, int score, Integer oldScore) {
		assert (oldScore == null || (int)oldScore != score); // we shouldn't be adding a node with the same score to the frontier twice
		// Check if the node is currently active with a lower score
		// if so, remove it from the frontier before we add the latest
		// score
		if (oldScore != null) {
			if (frontier.containsKey(oldScore)) {
				Set<SubgraphPathNode> oldLocation = frontier.get(oldScore);
				if (oldLocation.contains(node)) {
					oldLocation.remove(node);
					if (oldLocation.isEmpty()) {
						frontier.remove(oldScore);
					}
				}
			}
		}
		if (!frontier.containsKey(score)) {
			Set<SubgraphPathNode> value = Sets.newHashSet();
			frontier.put(score, value);
		}
		frontier.get(score).add(node);
	}
	private static SubgraphPathNode popFronter(NavigableMap<Integer, Set<SubgraphPathNode>> frontier) {
		Entry<Integer, Set<SubgraphPathNode>> entry = frontier.firstEntry();
		Set<SubgraphPathNode> value = entry.getValue();
		SubgraphPathNode head = value.iterator().next();
		value.remove(head);
		if (value.isEmpty()) {
			frontier.remove(entry.getKey());
		}
		return head;
	}
	public static final Function<SubgraphPathNode, Integer> SCORE_TOTAL_KMER = new Function<SubgraphPathNode, Integer>() {
		public Integer apply(SubgraphPathNode arg) {
			return arg.getWeight();
		}
	};
	public static final Function<SubgraphPathNode, Integer> SCORE_MAX_KMER = new Function<SubgraphPathNode, Integer>() {
		public Integer apply(SubgraphPathNode arg) {
			return arg.getMaxKmerWeight();
		}
	};
	private static Log expensiveSanityCheckLog = Log.getInstance(PathGraphAssembler.class);
	private boolean sanityCheckSubgraphs() {
		if (expensiveSanityCheckLog != null) {
			expensiveSanityCheckLog.warn("Expensive sanity checking is being performed. Performance will be poor");
			expensiveSanityCheckLog = null;
		}
		assert(subgraphs.size() == startingPaths.size());
		for (int i = 0; i < subgraphs.size(); i++) {
			boolean shouldBeReferenceStart = false;
			for (SubgraphPathNode n : subgraphs.get(i)) {
				for (SubgraphPathNode adj : prevPath(n)) {
					if (adj.containsNonReferenceKmer()) {
						assert(subgraphs.get(i).contains(adj));
					} else if (adj.containsReferenceKmer()) {
						shouldBeReferenceStart = true;
					}
				}
				for (SubgraphPathNode adj : nextPath(n)) {
					if (adj.containsNonReferenceKmer()) {
						assert(subgraphs.get(i).contains(adj));
					}
				}
				for (int j = 0; j < subgraphs.size(); j++) {
					if (i != j) {
						// disjoint
						assert(!subgraphs.get(j).contains(n));
					}
				}
			}
			if (shouldBeReferenceStart) {
				for (SubgraphPathNode n : startingPaths.get(i)) {
					boolean hasPrevReference = false;
					for (SubgraphPathNode adj : prevPath(n)) {
						if (adj.containsReferenceKmer()) {
							hasPrevReference = true;
							break;
						}
					}
					assert(hasPrevReference);
				}
			} else {
				for (SubgraphPathNode n : startingPaths.get(i)) {
					assert(n.containsNonReferenceKmer());
				}
			}
		}
		return true;
	}
	private boolean sanityCheckMemoization(
			Map<SubgraphPathNode, Integer> bestScore,
			Map<SubgraphPathNode, List<SubgraphPathNode>> bestPath,
			NavigableMap<Integer, Set<SubgraphPathNode>> frontier,
			Function<SubgraphPathNode, Integer> scoringFunction,
			boolean scoreReference, boolean scoreNonReference) {
		if (expensiveSanityCheckLog != null) {
			expensiveSanityCheckLog.warn("Expensive sanity checking is being performed. Performance will be poor");
			expensiveSanityCheckLog = null;
		}
		for (Entry<SubgraphPathNode, Integer> entry : bestScore.entrySet()) {
			SubgraphPathNode node = entry.getKey();
			assert(bestScore.containsKey(node));
			assert(bestPath.containsKey(node));
			assert(bestScore.get(node) == getScore(bestPath.get(node), scoringFunction, scoreReference, scoreNonReference));
		}
		for (Entry<SubgraphPathNode, List<SubgraphPathNode>> entry : bestPath.entrySet()) {
			SubgraphPathNode node = entry.getKey();
			List<SubgraphPathNode> path = entry.getValue();
			assert(path.size() == Sets.newHashSet(path).size()); // ensure simple path
			assert(bestScore.containsKey(node));
			assert(bestPath.containsKey(node));
			assert(bestScore.get(node) == getScore(bestPath.get(node), scoringFunction, scoreReference, scoreNonReference));
		}
		Set<SubgraphPathNode> frontierNodes = Sets.newHashSet();
		int frontierNodeCount = 0;
		for (Entry<Integer, Set<SubgraphPathNode>> entry : frontier.entrySet()) {
			int score = entry.getKey();
			for (SubgraphPathNode node : entry.getValue()) {
				assert(bestScore.containsKey(node));
				assert(bestPath.containsKey(node));
				assert(frontier.get(score).contains(node));
				assert(bestScore.get(node) >= getScore(bestPath.get(node), scoringFunction, scoreReference, scoreNonReference));
			}
			frontierNodeCount += entry.getValue().size();
			frontierNodes.addAll(entry.getValue());
		}
		assert(frontierNodeCount == frontierNodes.size()); // no duplicates
		return true;
	}
}
