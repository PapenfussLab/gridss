package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Stack;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphUtil;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNode;
import au.edu.wehi.idsv.graph.PathNodeFactory;
import au.edu.wehi.idsv.graph.WeightedSequenceGraphNodeUtil;
import au.edu.wehi.idsv.visualisation.DeBruijnPathGraphExporter;
import au.edu.wehi.idsv.visualisation.StaticDeBruijnPathGraphGexfExporter;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
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
public class PathGraphAssembler<T, PN extends DeBruijnPathNode<T>> extends DeBruijnPathGraph<T, PN> {
	private static final Log log = Log.getInstance(PathGraphAssembler.class);
	private final AssemblyParameters parameters;
	private final List<NonReferenceSubgraph> subgraphs = Lists.newArrayList();
	private static int timeoutGraphsWritten = 0;
	private boolean timeoutReached = false;
	private boolean anchorTimeoutReached = false;
	private boolean misassemblyDetected = false;
	private int totalNodesTraversed = 0;
	private class NonReferenceSubgraph {
		/**
		 * Node set of the non-reference subgraphs
		 */
		public final List<PN> nodes = Lists.newArrayList();
		/**
		 * Nodes with a connected reference kmer preceding the given path
		 */
		public final List<PN> referenceAnchoredPreceding = Lists.newArrayList();
		/**
		 * Nodes with a connected reference kmer following the given path
		 */
		public final List<PN> referenceAnchoredPost = Lists.newArrayList();
		public NonReferenceSubgraph(PN seed, BitSet processed) {
			assert(!isReference(seed));
			visitAll(seed, processed);
		}
		private void visitAll(PN seed, BitSet processed) {
			Stack<PN> frontier = new Stack<PN>();
			nodes.add(seed);
			processed.set(seed.getNodeId());
			frontier.add(seed);
			while (!frontier.isEmpty()) {
				PN node = frontier.pop();
				boolean added = false;
				for (PN adj : prev(node)) {
					boolean adjIsRef = isReference(adj);
					if (!added && adjIsRef) {
						referenceAnchoredPreceding.add(node);
						added = true;
					}
					if (!adjIsRef && !processed.get(adj.getNodeId())) {
						nodes.add(adj);
						processed.set(adj.getNodeId());
						frontier.add(adj);
					}
				}
				added = false;
				for (PN adj : next(node)) {
					boolean adjIsRef = isReference(adj);
					if (!added && adjIsRef) {
						referenceAnchoredPost.add(node);
						added = true;
					}
					if (!adjIsRef && !processed.get(adj.getNodeId())) {
						nodes.add(adj);
						processed.set(adj.getNodeId());
						frontier.add(adj);
					}
				}
			}
		}
		private List<PN> traverseNonReference(
				Collection<PN> startNodes,
				Collection<PN> endNodes,
				boolean traverseForward
				) {
			PathGraphTraverse<T, PN> trav = new PathGraphTraverse<T, PN>(
					PathGraphAssembler.this,
					parameters.subgraphMaxPathTraversalNodes,
					traverseForward,
					false,
					parameters.subgraphAssemblyTraversalMaximumBranchingFactor,
					Integer.MAX_VALUE);
			List<PN> best;
			if (endNodes != null) {
				best = trav.traverse(startNodes, endNodes); 
			} else {
				best = trav.traverse(startNodes);
			}
			totalNodesTraversed += trav.getNodesTraversed();
			timeoutReached |= trav.getNodesTraversed() >= parameters.subgraphMaxPathTraversalNodes;
			return best;
		}
		/**
		 * Generates a reference anchor contig
		 * @param node non-reference starting node
		 * @param traverseForward direction of traversal
		 * @param targetLength maximum number of non-reference kmers to traverse
		 * @return reference path
		 */
		private List<PN> traverseReference(PN node, boolean traverseForward, int targetLength) {
			assert(!isReference(node));
			PathGraphTraverse<T, PN> trav = new PathGraphTraverse<T, PN>(PathGraphAssembler.this, parameters.subgraphMaxPathTraversalNodes,
					traverseForward, true,
					// if anchor has previously timed out on a different breakend subgraph
					// just do a greedy traversal as we're likely to time out again
					anchorTimeoutReached ? 1 : parameters.subgraphAssemblyTraversalMaximumBranchingFactor,
					// assemble anchorAssemblyLength bases / at least the length of the breakend
					node.length() + targetLength);
			List<PN> anchor = trav.traverse(ImmutableList.of(node));
			totalNodesTraversed += trav.getNodesTraversed();
			anchorTimeoutReached |= trav.getNodesTraversed() >= parameters.subgraphMaxPathTraversalNodes; 
			timeoutReached |= anchorTimeoutReached;
			assert(anchor != null); // assembly includes starting node so must include that
			assert(anchor.size() > 0); // assembly includes starting node so must include that
			anchor.remove(traverseForward ? 0 : anchor.size() - 1); // don't repeat the non-reference node we used as the anchor
			return anchor;
		}
		/**
		 * Returns the highest scoring contig from the non-reference subgraph 
		 * @param scoringFunction node scoring function
		 * @return non-reference contig
		 */
		private List<PN> assembleSubgraph() {
			List<PN> best = null;
			if (!referenceAnchoredPreceding.isEmpty() && referenceAnchoredPreceding.isEmpty()) {
				if (referenceAnchoredPreceding.size() >= referenceAnchoredPreceding.size()) {
					best = traverseNonReference(referenceAnchoredPreceding, referenceAnchoredPost, true);
				} else {
					// going backwards is faster since we have fewer end paths than start paths
					best = traverseNonReference(referenceAnchoredPost, referenceAnchoredPreceding, false);
				}
				if (best == null) {
					//    A - B - C
					//     \     /
					//  R - R - R - R
					//
					// No path from C to A results in null contig
					// TODO: assemble both AR and CR as two separate contigs instead of just
					// falling through to an RC assembly
					misassemblyDetected = true;
					log.debug(String.format("Breakpoint assembly failed at %s. Falling back to breakend assembly.",
							referenceAnchoredPreceding.get(0).toString()));
				}
			}
			if (best == null && !referenceAnchoredPreceding.isEmpty()) {
				best = traverseNonReference(referenceAnchoredPreceding, null, true);
			}
			if (best == null && !referenceAnchoredPost.isEmpty()) {
				best = traverseNonReference(referenceAnchoredPost, null, false);
			}
			if (best == null) {
				// no anchors/anchored assemby failed -> just assemble what we can
				best = traverseNonReference(nodes, null, true);
			}
			assert(best != null); // subgraphs contain at least one node
			assert(best.size() > 0);
			
			int targetAnchorAssemblyLength = Math.max(parameters.anchorAssemblyLength, WeightedSequenceGraphNodeUtil.nodeLength(best));
			List<PN> preAnchor = traverseReference(best.get(0), false, targetAnchorAssemblyLength);			
			List<PN> postAnchor = traverseReference(best.get(best.size() - 1), true, targetAnchorAssemblyLength);
			ArrayList<PN> result = new ArrayList<PN>(best.size() + preAnchor.size() + postAnchor.size());
			result.addAll(preAnchor);
			result.addAll(best);
			result.addAll(postAnchor);
			assertHasSingleNonReferenceKmerChain_PN(result);
			return result;
		}
	}
	public PathGraphAssembler(DeBruijnGraph<T> graph, AssemblyParameters parameters, Collection<T> seeds, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		super(graph, seeds, factory, tracker);
		this.parameters = parameters;
	}
	public List<List<PN>> assembleContigs() {
		return assembleContigs(null);
	}
	public List<List<PN>> assembleContigs(DeBruijnPathGraphExporter<T, PN> graphExporter) {
		splitOutReferencePaths();
		calcNonReferenceSubgraphs();
		return assemblyNonReferenceContigs(graphExporter);
	}
	private boolean assertHasSingleNonReferenceKmerChain_PN(List<PN> contig) {
		int chainStarts = isReference(contig.get(0)) ? 0 : 1;
		for (int i = 1; i < contig.size(); i++) {
			if (isReference(contig.get(i - 1)) && !isReference(contig.get(i))) {
				chainStarts++;
			}
		}
		assert(chainStarts == 1);
		return chainStarts == 1;
	}
	private List<List<PN>> assemblyNonReferenceContigs(DeBruijnPathGraphExporter<T, PN> graphExporter) {
		List<List<PN>> result = Lists.newArrayList();
		for (NonReferenceSubgraph nrsg: subgraphs) {
			List<PN> contig = nrsg.assembleSubgraph();
			result.add(contig);
		}
		// return the best assembly first
		final DeBruijnGraph<PN> g = this;
		Collections.sort(result, new Ordering<List<PN>>() {
			@Override
			public int compare(List<PN> arg0, List<PN> arg1) {
				return Ints.compare(
						DeBruijnGraphUtil.scorePath(g, arg1, false, true),
						DeBruijnGraphUtil.scorePath(g, arg0, false, true));
			}}); 
		if (graphExporter != null) {
			graphExporter.snapshot(this);
			//graphExporter.annotateStartingPaths(startingPaths);
			graphExporter.contigs(result);
		} else if ((misassemblyDetected || timeoutReached) && parameters.debruijnGraphVisualisationDirectory != null) {
			// weren't planning to export but we encountering something worth logging
			parameters.debruijnGraphVisualisationDirectory.mkdirs();
			graphExporter = new StaticDeBruijnPathGraphGexfExporter<T, PN>();
			graphExporter.snapshot(this);
			//graphExporter.annotateStartingPaths(startingPaths);
			graphExporter.contigs(result);
			graphExporter.saveTo(new File(parameters.debruijnGraphVisualisationDirectory,
					String.format("pathTraversalTimeout-%d.subgraph.gexf", ++timeoutGraphsWritten)));
		}
		tracker.assemblyNonReferenceContigs(result, totalNodesTraversed);
		return result;
	}
	/**
	 * Separates non-reference kmers of the subgraph into disjoint non-reference subgraphs 
	 */
	private void calcNonReferenceSubgraphs() {
		BitSet processed = new BitSet(getMaxNodeId() + 1);
		for (PN node : getPaths()) {
			if (!processed.get(node.getNodeId()) && !isReference(node)) {
				// start of a new non-reference subgraph
				NonReferenceSubgraph nrsg = new NonReferenceSubgraph(node, processed); 
				subgraphs.add(nrsg);
			}
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheckSubgraphs());
		}
		//tracker.calcNonReferenceSubgraphs(subgraphs, startingPaths);
	}
	private static Log expensiveSanityCheckLog = Log.getInstance(PathGraphAssembler.class);
	private boolean sanityCheckSubgraphs() {
		if (expensiveSanityCheckLog != null) {
			expensiveSanityCheckLog.warn("Expensive sanity checking is being performed. Performance will be poor");
			expensiveSanityCheckLog = null;
		}
		for (NonReferenceSubgraph subgraph : subgraphs) {
			for (PN n : subgraph.nodes) {
				assert(subgraph.nodes.containsAll(subgraph.referenceAnchoredPreceding));
				assert(subgraph.nodes.containsAll(subgraph.referenceAnchoredPost));
				assert(!isReference(n));
				for (PN adj : prev(n)) {
					if (isReference(adj)) {
						assert(subgraph.referenceAnchoredPreceding.contains(n));
					} else {
						subgraph.nodes.contains(adj);
					}
				}
				for (PN adj : next(n)) {
					if (isReference(adj)) {
						assert(subgraph.referenceAnchoredPost.contains(n));
					} else {
						subgraph.nodes.contains(adj);
					}
				}
			}
		}
		List<PN> allSubgraphs = Lists.newArrayList(Iterables.concat(Iterables.transform(subgraphs, new Function<NonReferenceSubgraph, Iterable<PN>>() {
			@Override
			public Iterable<PN> apply(NonReferenceSubgraph input) {
				return input.nodes;
			}
		})));
		HashSet<PN> allSubgraphsNodeSet = Sets.newHashSet(allSubgraphs);
		// PNs occurs exactly once
		assert(allSubgraphs.size() == allSubgraphsNodeSet.size());
		// all non-reference PNs are included in subgraphs
		assert(allSubgraphsNodeSet.size() == Iterables.size(Iterables.filter(allNodes(), new Predicate<PN>() {
			@Override
			public boolean apply(PN input) {
				return !isReference(input);
			}
		})));
		return true;
	}
}
