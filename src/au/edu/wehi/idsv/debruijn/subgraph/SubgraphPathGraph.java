package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.LinkedList;

import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;

public class SubgraphPathGraph extends DeBruijnPathGraph<DeBruijnSubgraphNode, SubgraphPathNode> {
	public SubgraphPathGraph(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, long seed) {
		super(graph, seed);
	}
	@Override
	protected SubgraphPathNode createPathNode(LinkedList<Long> path, DeBruijnGraphBase<DeBruijnSubgraphNode> graph) {
		return new SubgraphPathNode(path, graph);
	}
}