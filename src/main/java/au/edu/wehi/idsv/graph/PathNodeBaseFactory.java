package au.edu.wehi.idsv.graph;

import java.util.Collection;

import com.google.common.collect.ImmutableList;

public class PathNodeBaseFactory<T> implements PathNodeFactory<T, PathNode<T>> {
	public final WeightedDirectedGraph<T> graph;
	public PathNodeBaseFactory(WeightedDirectedGraph<T> graph) {
		this.graph = graph;
	}
	@Override
	public PathNode<T> splitPathNode(PathNode<T> unsplit, int startIndex, int length) {
		return new PathNode<T>(ImmutableList.of(unsplit), startIndex, length, graph);
	}
	@Override
	public PathNode<T> concatPathNodes(Iterable<PathNode<T>> path) {
		return PathNode.concat(path, graph);
	}
	@Override
	public PathNode<T> createPathNode(Collection<T> path) {
		return new PathNode<T>(path, graph);
	}
}
