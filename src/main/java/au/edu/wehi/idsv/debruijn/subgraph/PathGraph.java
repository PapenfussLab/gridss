package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

public class PathGraph extends DeBruijnPathGraph<DeBruijnSubgraphNode, SubgraphPathNode> {
	public PathGraph(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, long seed, SubgraphAssemblyAlgorithmTracker tracker) {
		super(graph, seed, new SubgraphPathNodeFactory(), tracker);
	}
	private int referencePathsSplit = 0;
	public int getReferencePathsSplit() { return referencePathsSplit; }
	/**
	 * Breaks paths into paths in which all kmers are either reference or non-reference kmers 
	 */
	public void splitOutReferencePaths() {
		// get full list of splits required upfront (since we'll be changing the node set as we iterate) 
		List<SubgraphPathNode> toSplit = Lists.newArrayList(Iterables.filter(getPaths(), new Predicate<SubgraphPathNode>() {
			public boolean apply(SubgraphPathNode arg) {
				return arg.containsReferenceKmer() && arg.containsNonReferenceKmer();
			}
		}));
		for (SubgraphPathNode n : toSplit) {
			List<Integer> lengths = new ArrayList<>();
			List<Long> path = n.getPath();
			int currentLength = 0;
			boolean currentIsReference = getGraph().getKmer(path.get(0)).isReference();
			for (long kmer : path) {
				boolean isRef = getGraph().getKmer(kmer).isReference();
				if (currentIsReference == isRef) {
					currentLength++;
				} else {
					currentIsReference = isRef;
					lengths.add(currentLength);
					currentLength = 1;
				}
			}
			lengths.add(currentLength);
			split(n, lengths);
			referencePathsSplit += lengths.size() - 1;
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(assertReferenceKmersSplit());
		}
	}
	private boolean assertReferenceKmersSplit() {
		for (SubgraphPathNode n : getPaths()) {
			assert(n.containsNonReferenceKmer() != n.containsReferenceKmer());
			for (long kmer : n.getPath()) {
				assert(getGraph().getKmer(kmer).isReference() == n.containsReferenceKmer());
			}
		}
		return true;
	}
}
