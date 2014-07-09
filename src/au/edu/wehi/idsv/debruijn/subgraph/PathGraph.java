package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.List;

import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

public class PathGraph extends DeBruijnPathGraph<DeBruijnSubgraphNode, SubgraphPathNode> {
	public PathGraph(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, long seed) {
		super(graph, seed, new SubgraphPathNodeFactory());
	}
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
			List<Integer> lengths = Lists.newArrayList();
			List<Long> path = n.getPath();
			int currentLength = 0;
			boolean currentIsReference = graph.getKmer(path.get(0)).isReference();
			for (long kmer : path) {
				boolean isRef = graph.getKmer(kmer).isReference();
				if (currentIsReference == isRef) {
					currentLength++;
				} else {
					currentIsReference = isRef;
					lengths.add(currentLength);
					currentLength = 0;
				}
				lengths.add(currentLength);
			}
			split(n, lengths);
		}
	}
}
