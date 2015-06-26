package au.edu.wehi.idsv.debruijn;

import java.util.Collection;

import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.graph.WeightedDirectedGraph;

public class DeBruijnPathNode<T> extends PathNode<T> implements DeBruijnSequenceGraphNode {
	private int refCount = 0;
	private void addRefCounts(Iterable<T> nodes) {
		refCount = 0;
		for (T n : nodes) {
			addRefCounts(n);
		}
	}
	private void addRefCounts(T node) {
		if (getGraph().isReference(node)) {
			refCount++;
		}
	}
	public DeBruijnPathNode(Collection<T> nodes, DeBruijnGraph<T> graph) {
		super(nodes, graph);
		addRefCounts(nodes);
	}
	public DeBruijnPathNode(Iterable<? extends DeBruijnPathNode<T>> nodes, int startOffset, int length, DeBruijnGraph<T> graph) {
		super(nodes, startOffset, length, graph);
		addRefCounts(allNodes());
	}
	@Override
	protected void merge(int offset, PathNode<T> pn, int pnOffset, WeightedDirectedGraph<T> graph) {
		super.merge(offset, pn, pnOffset, graph);
		addRefCounts(pn.getPathAllNodes().get(pnOffset));
	}
	public boolean isReference() {
		return refCount > 0;
	}
	private DeBruijnGraph<T> getGraph() { return (DeBruijnGraph<T>)super.graph; }
	@Override
	public long kmer(int offset) {
		return getGraph().getKmer(this.getPath().get(offset));
	}
}
