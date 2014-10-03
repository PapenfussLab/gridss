package au.edu.wehi.idsv.debruijn.subgraph;

import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.PathNodeFactory;

/**
 * Factory for generating SubgraphPathNode instances
 * @author Daniel Cameron
 *
 */
public class SubgraphPathNodeFactory implements PathNodeFactory<DeBruijnSubgraphNode, SubgraphPathNode> {
	public SubgraphPathNode createPathNode(Iterable<Long> path, DeBruijnGraphBase<DeBruijnSubgraphNode> graph) {
		return new SubgraphPathNode(path, graph);
	}
	public SubgraphPathNode createPathNode(SubgraphPathNode unsplit, int startIndex, int length, DeBruijnGraphBase<DeBruijnSubgraphNode> graph) {
		return new SubgraphPathNode(unsplit, startIndex, length, graph);
	}
	public SubgraphPathNode createPathNode(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, Iterable<SubgraphPathNode> path) {
		return new SubgraphPathNode(graph, path);
	}
}