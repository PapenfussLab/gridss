package au.edu.wehi.idsv.graph;

import java.util.Collection;

public class PathNodeBaseFactory<T> implements PathNodeFactory<T, PathNode<T>> {
	@Override
	public PathNode<T> splitPathNode(PathNode<T> unsplit, int startIndex, int length, WeightedDirectedGraph<T> graph) {
		return new PathNode<T>(unsplit, startIndex, length, graph);
	}
	@Override
	public PathNode<T> concatPathNodes(Iterable<PathNode<T>> path, WeightedDirectedGraph<T> graph) {
		return PathNode.concat(path, graph);
	}
	@Override
	public PathNode<T> createPathNode(Collection<T> path, WeightedDirectedGraph<T> graph) {
		return new PathNode<T>(path, graph);
	}
}
