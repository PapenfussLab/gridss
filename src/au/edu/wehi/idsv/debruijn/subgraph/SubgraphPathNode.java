package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.LinkedList;

import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.PathNode;

public class SubgraphPathNode extends PathNode<DeBruijnSubgraphNode> {
	public SubgraphPathNode(LinkedList<Long> path, DeBruijnGraphBase<DeBruijnSubgraphNode> graph) {
		super(path, graph);
	}
}
