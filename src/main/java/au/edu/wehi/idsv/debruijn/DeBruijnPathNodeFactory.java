package au.edu.wehi.idsv.debruijn;

import java.util.Collection;

import au.edu.wehi.idsv.graph.PathNodeFactory;
import au.edu.wehi.idsv.graph.WeightedSequenceGraphNodeUtil;

import com.google.common.collect.ImmutableList;

public class DeBruijnPathNodeFactory<T> implements PathNodeFactory<T, DeBruijnPathNode<T>> {
	public final DeBruijnGraph<T> graph;
	public DeBruijnPathNodeFactory(DeBruijnGraph<T> graph) {
		this.graph = graph;
	}
	@Override
	public DeBruijnPathNode<T> splitPathNode(DeBruijnPathNode<T> unsplit, int startIndex, int length) {
		return new DeBruijnPathNode<T>(ImmutableList.of(unsplit), startIndex, length, graph);
	}
	@Override
	public DeBruijnPathNode<T> concatPathNodes(Iterable<DeBruijnPathNode<T>> path) {
		return new DeBruijnPathNode<T>(path, 0, WeightedSequenceGraphNodeUtil.nodeLength(path), graph);
	}
	@Override
	public DeBruijnPathNode<T> createPathNode(Collection<T> path) {
		return new DeBruijnPathNode<T>(path, graph);
	}
}
