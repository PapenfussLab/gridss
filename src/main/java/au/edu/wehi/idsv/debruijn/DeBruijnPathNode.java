package au.edu.wehi.idsv.debruijn;

import java.util.Collection;

import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.graph.WeightedDirectedGraph;

public class DeBruijnPathNode<T> extends PathNode<T> {
	private int refCount = 0;
	private void addRefCounts(Iterable<T> nodes, DeBruijnGraph<T> graph) {
		refCount = 0;
		for (T n : nodes) {
			addRefCounts(n, graph);
		}
	}
	private void addRefCounts(T node, DeBruijnGraph<T> graph) {
		if (graph.isReference(node)) {
			refCount++;
		}
	}
	public DeBruijnPathNode(Collection<T> nodes, DeBruijnGraph<T> graph) {
		super(nodes, graph);
		addRefCounts(nodes, graph);
	}
	public DeBruijnPathNode(Iterable<? extends DeBruijnPathNode<T>> nodes, int startOffset, int length, DeBruijnGraph<T> graph) {
		super(nodes, startOffset, length, graph);
		addRefCounts(getPath(), graph);
	}
	@Override
	protected void merge(int offset, PathNode<T> pn, int pnOffset, WeightedDirectedGraph<T> graph) {
		super.merge(offset, pn, pnOffset, graph);
		addRefCounts(pn.getPathAllNodes().get(pnOffset), (DeBruijnGraph<T>)graph);
	}
	public boolean isReference() {
		return refCount > 0;
	}
}
