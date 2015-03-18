package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.PathNode;
import au.edu.wehi.idsv.visualisation.DeBruijnPathGraphExporter;
import au.edu.wehi.idsv.visualisation.StaticDeBruijnSubgraphPathGraphGexfExporter;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

/**
 * Reduces the de bruijn graph to a set of paths before assembling
 * @author Daniel Cameron
 *
 */
public class PathGraphAssembler extends PathGraph {
	//private static final Log log = Log.getInstance(PathGraphAssembler.class);
	private final AssemblyParameters parameters;
	//private final List<Set<SubgraphPathNode>> subgraphs = Lists.newArrayList();
	private final List<Iterable<SubgraphPathNode>> startingPaths = Lists.newArrayList();
	private static int timeoutGraphsWritten = 0;
	private boolean timeoutReached = false;
	private int totalNodesTraversed = 0;
	public PathGraphAssembler(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, AssemblyParameters parameters, long seed, SubgraphAssemblyAlgorithmTracker tracker) {
		super(graph, seed, tracker);
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
		final Function<SubgraphPathNode, Integer> scoringFunction = SubgraphHelper.SCORE_TOTAL_KMER;
		List<List<SubgraphPathNode>> result = Lists.newArrayList();
		for (Iterable<SubgraphPathNode> starts : startingPaths) {
			List<SubgraphPathNode> contig = assembleSubgraph(starts, scoringFunction);
			result.add(contig);
		}
		// return the best assembly first
		Collections.sort(result, new Ordering<List<SubgraphPathNode>>() {
			@Override
			public int compare(List<SubgraphPathNode> arg0, List<SubgraphPathNode> arg1) {
				return Ints.compare(
						SubgraphHelper.getScore(arg1, scoringFunction, false, true),
						SubgraphHelper.getScore(arg0, scoringFunction, false, true));
			}}); 
		if (graphExporter != null) {
			graphExporter.snapshot(this);
			graphExporter.annotateStartingPaths(startingPaths);
			graphExporter.contigs(result);
		} else if (timeoutReached && parameters.debruijnGraphVisualisationDirectory != null) {
			// weren't planning to export but we timed out
			parameters.debruijnGraphVisualisationDirectory.mkdirs();
			graphExporter = new StaticDeBruijnSubgraphPathGraphGexfExporter(this.parameters.k);
			graphExporter.snapshot(this);
			graphExporter.annotateStartingPaths(startingPaths);
			graphExporter.contigs(result);
			graphExporter.saveTo(new File(parameters.debruijnGraphVisualisationDirectory,
					String.format("pathTraversalTimeout-%d.subgraph.gexf", ++timeoutGraphsWritten)));
		}
		List<LinkedList<Long>> resultKmers = Lists.newArrayList(Iterables.transform(result, new Function<List<SubgraphPathNode>, LinkedList<Long>>() {
			@Override
			public LinkedList<Long> apply(List<SubgraphPathNode> arg0) {
				return Lists.newLinkedList(Lists.newArrayList(SubgraphPathNode.kmerIterator(arg0)));
			}
		}));
		tracker.assemblyNonReferenceContigs(result, resultKmers, totalNodesTraversed);
		return resultKmers;
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
				toProcess.removeAll(newNonreferenceSubgraph(node));
			}
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheckSubgraphs());
		}
		//tracker.calcNonReferenceSubgraphs(subgraphs, startingPaths);
	}
	private Set<SubgraphPathNode> newNonreferenceSubgraph(SubgraphPathNode seed) {
		Set<SubgraphPathNode> nodes = Sets.newHashSet();
		List<SubgraphPathNode> startingNodes = new ArrayList<SubgraphPathNode>();
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
			startingNodes = Lists.newArrayListWithCapacity(nodes.size());
			startingNodes.addAll(nodes);
		}
		//subgraphs.add(nodes);
		startingPaths.add(startingNodes);
		return nodes;
	}
	/**
	 * Returns the best possible non-reference contig starting from one of the starting nodes. 
	 * @param startingNodes available starting nodes
	 * @param scoringFunction node scoring function
	 * @return non-reference contig, preceeded by reference assembly 
	 */
	private List<SubgraphPathNode> assembleSubgraph(
			Iterable<SubgraphPathNode> startingNodes,
			Function<SubgraphPathNode, Integer> scoringFunction) {
		assert(startingNodes != null);
		int branchingFactor = parameters.subgraphAssemblyTraversalMaximumBranchingFactor;
		
		PathGraphTraverse trav = new PathGraphTraverse(this, parameters.maxPathTraversalNodes, scoringFunction, true, false, branchingFactor, Integer.MAX_VALUE);
		List<SubgraphPathNode> best = trav.traverse(startingNodes);
		totalNodesTraversed += trav.getNodesTraversed();
		timeoutReached |= trav.getNodesTraversed() >= parameters.maxPathTraversalNodes;
		assert(best != null); // subgraphs contain at least one node
		assert(best.size() > 0);
		SubgraphPathNode breakendAnchorStart = best.get(0);
		PathGraphTraverse anchorTrav = new PathGraphTraverse(this, parameters.maxPathTraversalNodes, scoringFunction, false, true, branchingFactor,
				breakendAnchorStart.length() + Math.max(parameters.anchorAssemblyLength, PathNode.kmerLength(best)));
		List<SubgraphPathNode> anchor = anchorTrav.traverse(ImmutableList.of(breakendAnchorStart));
		totalNodesTraversed += trav.getNodesTraversed();
		timeoutReached |= trav.getNodesTraversed() >= parameters.maxPathTraversalNodes;
		assert(anchor != null); // assembly includes starting node so must include that
		assert(anchor.size() > 0); // assembly includes starting node so must include that
		anchor.remove(anchor.size() - 1); // don't repeat the non-reference node we used as the anchor
		ArrayList<SubgraphPathNode> result = new ArrayList<SubgraphPathNode>(best.size() + anchor.size());
		result.addAll(anchor);
		result.addAll(best);
		return result;
	}
	private static Log expensiveSanityCheckLog = Log.getInstance(PathGraphAssembler.class);
	private boolean sanityCheckSubgraphs() {
		if (expensiveSanityCheckLog != null) {
			expensiveSanityCheckLog.warn("Expensive sanity checking is being performed. Performance will be poor");
			expensiveSanityCheckLog = null;
		}
		/*
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
		*/
		return true;
	}
}
