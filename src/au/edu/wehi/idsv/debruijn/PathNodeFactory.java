package au.edu.wehi.idsv.debruijn;

public interface PathNodeFactory<T extends DeBruijnNodeBase, PN extends PathNode<T>> {
	PN createPathNode(Iterable<Long> path, DeBruijnGraphBase<T> graph);
	PN createPathNode(PN unsplit, int startIndex, int length, DeBruijnGraphBase<T> graph);
	PN createPathNode(DeBruijnGraphBase<T> graph, Iterable<PN> path);
}
