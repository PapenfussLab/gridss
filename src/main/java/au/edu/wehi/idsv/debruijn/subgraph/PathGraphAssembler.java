package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphUtil;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.graph.PathNodeFactory;
import au.edu.wehi.idsv.visualisation.DeBruijnPathGraphExporter;
import au.edu.wehi.idsv.visualisation.StaticDeBruijnPathGraphGexfExporter;
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
public class PathGraphAssembler<T, PN extends PathNode<T>> extends DeBruijnPathGraph<T, PN> {
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
		public final Set<PN> nodes = Sets.newHashSet();
		/**
		 * Nodes with a connected reference kmer preceding the given path
		 */
		public final List<PN> referenceAnchoredPreceding = Lists.newArrayList();
		/**
		 * Nodes with a connected reference kmer following the given path
		 */
		public final List<PN> referenceAnchoredPost = Lists.newArrayList();
		public NonReferenceSubgraph(PN seed) {
			assert(!isReference(seed));
			visitAll(seed);
		}
		private void visitAll(PN seed) {
			Queue<PN> frontier = new ArrayDeque<PN>();
			nodes.add(seed);
			frontier.add(seed);
			while (!frontier.isEmpty()) {
				PN node = frontier.poll();
				boolean added = false;
				for (PN prev : prev(node)) {
					boolean adjIsRef = isReference(prev);
					if (!added && adjIsRef) {
						referenceAnchoredPreceding.add(node);
						added = true;
					}
					if (!adjIsRef && !nodes.contains(prev)) {
						nodes.add(prev);
						frontier.add(prev);
					}
				}
				added = false;
				for (PN next : next(node)) {
					boolean adjIsRef = isReference(next);
					if (!added && adjIsRef) {
						referenceAnchoredPost.add(node);
						added = true;
					}
					if (!adjIsRef && !nodes.contains(next)) {
						nodes.add(next);
						frontier.add(next);
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
					parameters.maxPathTraversalNodes,
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
			timeoutReached |= trav.getNodesTraversed() >= parameters.maxPathTraversalNodes;
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
			PathGraphTraverse<T, PN> trav = new PathGraphTraverse<T, PN>(PathGraphAssembler.this, parameters.maxPathTraversalNodes,
					traverseForward, true,
					// if anchor has previously timed out on a different breakend subgraph
					// just do a greedy traversal as we're likely to time out again
					anchorTimeoutReached ? 1 : parameters.subgraphAssemblyTraversalMaximumBranchingFactor,
					// assemble anchorAssemblyLength bases / at least the length of the breakend
					node.length() + targetLength);
			List<PN> anchor = trav.traverse(ImmutableList.of(node));
			totalNodesTraversed += trav.getNodesTraversed();
			anchorTimeoutReached |= trav.getNodesTraversed() >= parameters.maxPathTraversalNodes; 
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
			
			int targetAnchorAssemblyLength = Math.max(parameters.anchorAssemblyLength, PN.nodeLength(best));
			List<PN> preAnchor = traverseReference(best.get(0), false, targetAnchorAssemblyLength);			
			List<PN> postAnchor = traverseReference(best.get(best.size() - 1), true, targetAnchorAssemblyLength);
			ArrayList<PN> result = new ArrayList<PN>(best.size() + preAnchor.size() + postAnchor.size());
			result.addAll(preAnchor);
			result.addAll(best);
			result.addAll(postAnchor);
			return result;
		}
	}
	public PathGraphAssembler(DeBruijnGraph<T> graph, AssemblyParameters parameters, Collection<T> seeds, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		super(graph, seeds, factory, tracker);
		this.parameters = parameters;
	}
	public List<LinkedList<T>> assembleContigs() {
		return assembleContigs(null);
	}
	public List<LinkedList<T>> assembleContigs(DeBruijnPathGraphExporter<T, PN> graphExporter) {
		splitOutReferencePaths();
		calcNonReferenceSubgraphs();
		return assemblyNonReferenceContigs(graphExporter);
	}
	private List<LinkedList<T>> assemblyNonReferenceContigs(DeBruijnPathGraphExporter<T, PN> graphExporter) {
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
		List<LinkedList<T>> resultKmers = Lists.newArrayList(Iterables.transform(result, new Function<List<PN>, LinkedList<T>>() {
			@Override
			public LinkedList<T> apply(List<PN> arg0) {
				return Lists.newLinkedList(Lists.newArrayList(PN.nodeIterator(arg0.iterator())));
			}
		}));
		tracker.assemblyNonReferenceContigs(result, resultKmers, totalNodesTraversed);
		return resultKmers;
	}
	/**
	 * Separates non-reference kmers of the subgraph into disjoint non-reference subgraphs 
	 */
	private void calcNonReferenceSubgraphs() {
		Set<PN> toProcess = Sets.newHashSet(getPaths());
		while (!toProcess.isEmpty()) {
			PN node = toProcess.iterator().next();
			toProcess.remove(node);
			if (!isReference(node)) {
				NonReferenceSubgraph nrsg = new NonReferenceSubgraph(node); 
				subgraphs.add(nrsg);
				toProcess.removeAll(nrsg.nodes);
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
		/*
		assert(subgraphs.size() == startingPaths.size());
		for (int i = 0; i < subgraphs.size(); i++) {
			boolean shouldBeReferenceStart = false;
			for (PN n : subgraphs.get(i)) {
				for (PN adj : prev(n)) {
					if (adj.containsNonReferenceKmer()) {
						assert(subgraphs.get(i).contains(adj));
					} else if (adj.containsReferenceKmer()) {
						shouldBeReferenceStart = true;
					}
				}
				for (PN adj : next(n)) {
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
				for (PN n : startingPaths.get(i)) {
					boolean hasPrevReference = false;
					for (PN adj : prev(n)) {
						if (adj.containsReferenceKmer()) {
							hasPrevReference = true;
							break;
						}
					}
					assert(hasPrevReference);
				}
			} else {
				for (PN n : startingPaths.get(i)) {
					assert(n.containsNonReferenceKmer());
				}
			}
		}
		*/
		return true;
	}
}
