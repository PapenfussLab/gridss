package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.SortedSet;

import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.debruijn.PathNode;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;
/**
 * Traverses a given subgraph and iteratively generates contigs from the
 * putative SVs in a greedy fashion.
 * 
 * Each non-reference kmer can be included in only a single assembly.
 *  
 * @author Daniel Cameron
 *
 */
public class GreedySubgraphContigAssembler extends SubgraphPathContigAssembler {
	public GreedySubgraphContigAssembler(DeBruijnReadGraph graph, long seedKmer) {
		super(graph, seedKmer);
	}
	/**
	 * Extracts all SV-supporting contigs from the de bruijn graph by greedy assembly from the start kmer
	 * @return kmer contig in subgraph
	 */
	public List<LinkedList<Long>> assembleContigs() {
		pathGraph.removeSelfIntersectingPaths();
		//pathGraph.collapseSimilarPaths(2, true);
		List<LinkedList<Long>> contigs = Lists.newArrayList();
		SubgraphPathNode seed;
		while ((seed = getBestSeed()) != null) {
			LinkedList<Long> contig = assembleContig(seed);
			if (contig != null && contig.size() > 0) {
				contigs.add(contig);
			}
		}
		return contigs;
	}
	/**
	 * Assembles a contig from the remaining kmers in the subgraph
	 * @param seed starting kmer
	 * @return kmer contig
	 */
	private LinkedList<Long> assembleContig(SubgraphPathNode seed) {
		LinkedList<SubgraphPathNode> contigKmers = pathGraph.greedyTraverse(seed, pathGraph.ByMaxKmerWeightDesc, pathGraph.ByMaxKmerWeightDesc);
		LinkedList<Long> path = Lists.newLinkedList();
		Iterators.addAll(path, PathNode.kmerIterator(contigKmers));
		trimPathToExpectedStructure(path, seed);
		removeNonReferencePathKmers(path);
		return path;
	}
	/**
	 * Gets the best kmer to start assembly from
	 * @return highest weighted non-reference kmer, null if no valid starting kmers can be found
	 */
	private SubgraphPathNode getBestSeed() {
		if (pathGraph.getPaths().isEmpty()) return null;
		return Collections.max(pathGraph.getPaths(), pathGraph.ByMaxKmerWeightDesc);
	}
}
