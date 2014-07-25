package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.ArrayDeque;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.AssemblyParameters.ContigAssemblyOrder;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
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
		return generateNonOverlapingContig(parameters.assemblyOrder);
		
	}	
	private List<LinkedList<Long>> generateNonOverlapingContigs(ContigAssemblyOrder algorithm) {
		List<LinkedList<Long>> results = Lists.newArrayList();
		// split graph into non-reference subgraphs
		List<Set<SubgraphPathNode>> nonReferenceSubgraphs = getNonReferenceSubgraphs();
		// TODO: need to get starting node list as well
		for (Set<SubgraphPathNode> nodeSet : nonReferenceSubgraphs) {
			// make assembly call
			List<SubgraphPathNode> contig = assembleSubgraph(nodeSet, algorithm);
			assert(contig != null);
			if (contig != null) {
				results.add(Lists.newLinkedList(Lists.newArrayList(SubgraphPathNode.kmerIterator(contig))));
			}
		}
		return results;
	}
	private List<SubgraphPathNode> assembleSubgraph(Set<SubgraphPathNode> nodeSet, ContigAssemblyOrder algorithm) {
		Function<SubgraphPathNode, Integer> scoringFunction = SCORE_MAX_KMER;
		List<SubgraphPathNode> best = null;
		int bestScore = Integer.MIN_VALUE;
		for (SubgraphPathNode startingNode : getStartingCandidates(nodeSet)) {
			List<SubgraphPathNode> current;
			switch (algorithm) {
				case GreedyMaxKmer:
					current = greedyAssembly(startingNode, nodeSet, scoringFunction);
					break;
				case OptimalMaxKmer:
					throw new IllegalArgumentException("OptimalMaxKmer Not Yet Implemented");
				default:
					throw new IllegalArgumentException("Unknown assembly order");
			}
			int currentScore = getScore(current, scoreMaxKmer);
			if (currentScore >= bestScore) {
				bestScore = currentScore;
				best = current;
				// traverse back to reference from starting node
			}
		}
		return best;
	}
	private int getScore(List<SubgraphPathNode> current, Function<SubgraphPathNode, Integer> scoringFunction) {
		int sum = 0;
		for (SubgraphPathNode n : current) {
			sum += scoringFunction.apply(current);
		}
		return sum;
	}
	private static final Function<SubgraphPathNode, Integer> SCORE_MAX_KMER = new Function<SubgraphPathNode, Integer>() {
		public Integer apply(SubgraphPathNode arg) {
			return arg.getMaxKmerWeight();
		}
	};
	private List<LinkedList<Long>> greedyAssembly(AssemblyParameters parameters) {
		PriorityQueue<SubgraphPathNode> seeds = new PriorityQueue<SubgraphPathNode>(16, pathGraph.ByMaxKmerWeightDesc);
		seeds.addAll(Sets.filter(pathGraph.getPaths(), new Predicate<SubgraphPathNode>() {
			public boolean apply(SubgraphPathNode arg) {
				boolean isNonReference = !arg.containsReferenceKmer(); 
				return isNonReference;
			}
		}));
		List<LinkedList<Long>> contigs = Lists.newArrayList();
		Set<SubgraphPathNode> visited = Sets.newHashSet();
		while (contigs.size() < parameters.maxContigsPerAssembly) {
			// flush out starting positions we've already visited
			while (!seeds.isEmpty() && visited.contains(seeds.peek())) seeds.poll();
			if (seeds.isEmpty()) break; // no more starting position
			// Found somewhere to start - off we go
			List<SubgraphPathNode> contig = greedyAssembly(seeds.poll(), visited);
			addNonReferenceReachable(visited, getFirstNonReference(contig));
			//addNonReferenceOnPath(visited, contig);
			contigs.add(Lists.newLinkedList(Lists.newArrayList(SubgraphPathNode.kmerIterator(contig))));
		}
		return contigs;
	}
	private void addNonReferenceOnPath(Set<SubgraphPathNode> visited, List<SubgraphPathNode> contig) {
		for (SubgraphPathNode node : contig) {
			if (node.containsNonReferenceKmer()) {
				assert(!node.containsReferenceKmer());
				// don't revisit non-reference kmers: anchoring two contigs to
				// nearby reference locations is ok
				visited.add(node);
			}
		}
	}
	private SubgraphPathNode getFirstNonReference(List<SubgraphPathNode> contig) {
		for (SubgraphPathNode node : contig) {
			if (node.containsNonReferenceKmer()) {
				assert(!node.containsReferenceKmer());
				return node;
			}
		}
		return null;
	}
	/**
	 * Adds all non-reference paths reachable from the given starting node to the visited set
	 * <p>Note: does not traverse elements already in the visited set </p> 
	 * @param visited set of visited nodes
	 * @param node node to start traversal from
	 */
	private void addNonReferenceReachable(Set<SubgraphPathNode> visited, SubgraphPathNode node) {
		if (node == null) return;
		ArrayDeque<SubgraphPathNode> frontier = new ArrayDeque<>();
		frontier.push(node);
		visited.add(node);
		while (!frontier.isEmpty()) {
			SubgraphPathNode n = frontier.poll();
			for (SubgraphPathNode adj : Iterables.concat(pathGraph.prevPath(n), pathGraph.nextPath(n))) {
				if (adj.containsReferenceKmer()) continue;
				if (visited.contains(adj)) continue;
				frontier.add(adj);
				visited.add(adj);
			}
		}
	}
	/**
	 * Assembles a contig from the remaining kmers in the subgraph
	 * @param seed starting kmer
	 * @param visited 
	 * @return kmer contig
	 */
	public LinkedList<SubgraphPathNode> greedyAssembly(SubgraphPathNode seed, Set<SubgraphPathNode> visited, Ordering<SubgraphPathNode> choice) {
		LinkedList<SubgraphPathNode> contigKmers = pathGraph.greedyTraverse(seed, choice, choice, visited);
		trimPath(contigKmers, seed);
		return contigKmers;
	}
	public LinkedList<SubgraphPathNode> memoizedAssembly(SubgraphPathNode seed, Set<SubgraphPathNode> visited) {
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
