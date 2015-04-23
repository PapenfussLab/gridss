package au.edu.wehi.idsv.graph;

import java.util.Collection;

public interface PathNodeFactory<T, PN extends PathNode<T>> {
	PN createPathNode(Collection<T> path, WeightedDirectedGraph<T> graph);
	PN splitPathNode(PN unsplit, int startIndex, int length, WeightedDirectedGraph<T> graph);
	PN concatPathNodes(Iterable<PN> path, WeightedDirectedGraph<T> graph);
}
