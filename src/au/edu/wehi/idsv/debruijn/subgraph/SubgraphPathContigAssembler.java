package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;
import java.util.SortedSet;

import au.edu.wehi.idsv.AssemblyParameters;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Reduces the de bruijn graph to a set of paths before assembling
 * @author Daniel Cameron
 *
 */
public class SubgraphPathContigAssembler {
	protected final DeBruijnReadGraph graph;
	protected final PathGraph pathGraph;
	protected final long seedKmer;
	public SubgraphPathContigAssembler(DeBruijnReadGraph graph, long seedKmer) {
		this.graph = graph;
		this.seedKmer = seedKmer;
		this.pathGraph = new PathGraph(graph, seedKmer);
	}
	public List<LinkedList<Long>> assembleContigs(AssemblyParameters parameters) {
		if (parameters.maxBaseMismatchForCollapse > 0) {
			pathGraph.collapseSimilarPaths(parameters.maxBaseMismatchForCollapse, parameters.collapseBubblesOnly);
		}
		pathGraph.splitOutReferencePaths();
		switch (parameters.assemblyOrder) {
			case GreedyMaxKmer:
				return greedyAssembly(parameters);
			case OptimalMaxKmer:
				throw new IllegalArgumentException("OptimalMaxKmer Not Yet Implemented");
			default:
				throw new IllegalArgumentException("Unknown assembly order");
		}
	}	
	private List<LinkedList<Long>> greedyAssembly(AssemblyParameters parameters) {
		SortedSet<SubgraphPathNode> seeds = Sets.newTreeSet(pathGraph.ByMaxKmerWeightDesc);
		// why doesn't Set have addAll(Iterable)? :(
		seeds.addAll(Lists.newArrayList(Iterables.filter(pathGraph.getPaths(), new Predicate<SubgraphPathNode>() {
			public boolean apply(SubgraphPathNode arg) {
				return !arg.containsReferenceKmer();
			}
		})));
		List<LinkedList<Long>> contigs = Lists.newArrayList();
		Set<SubgraphPathNode> visited = Sets.newHashSet();
		while (!seeds.isEmpty() && contigs.size() > parameters.maxContigsPerAssembly) {
			List<SubgraphPathNode> contig = greedyAssembly(seeds.first(), visited);
			for (SubgraphPathNode node : contig) {
				if (seeds.contains(node)) {
					visited.add(node);
				}
			}
			seeds.removeAll(contig);
			contigs.add(Lists.newLinkedList(Lists.newArrayList(SubgraphPathNode.kmerIterator(contig))));
		}
		return contigs;
	}
	/**
	 * Assembles a contig from the remaining kmers in the subgraph
	 * @param seed starting kmer
	 * @param visited 
	 * @return kmer contig
	 */
	public LinkedList<SubgraphPathNode> greedyAssembly(SubgraphPathNode seed, Set<SubgraphPathNode> visited) {
		LinkedList<SubgraphPathNode> contigKmers = pathGraph.greedyTraverse(seed, pathGraph.ByMaxKmerWeightDesc, pathGraph.ByMaxKmerWeightDesc, visited);
		trimPath(contigKmers, seed);
		return contigKmers;
	}
	/**
	 * Trims the path such that once the path is anchored to the reference,
	 * the path does not deviated from the reference   
	 * @param contigKmers
	 */
	public void trimPath(LinkedList<SubgraphPathNode> contigKmers, SubgraphPathNode seed) {
		ListIterator<SubgraphPathNode> it = contigKmers.listIterator();
		// advance forward to seed
		while (it.hasNext() && it.next() != seed);
		// go back to the reference anchor and consume the entire anchor
		boolean inReference = false;
		while (it.hasPrevious()) {
			SubgraphPathNode node = it.previous();
			if (node.containsReferenceKmer()) {
				inReference = true;
			} else if (inReference && !node.containsReferenceKmer()) {
				it.next();
				break;
			}
		}
		if (inReference) {
			// if we indeed had an reference anchor, drop everything before it
			while (it.hasPrevious()) {
				it.previous();
				it.remove();
			}
		}
	}
}
