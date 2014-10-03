package au.edu.wehi.idsv.debruijn;

public class PathNodeBaseFactory implements PathNodeFactory<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>> {
	@Override
	public PathNode<DeBruijnNodeBase> createPathNode(Iterable<Long> path, DeBruijnGraphBase<DeBruijnNodeBase> graph) {
		return new PathNode<DeBruijnNodeBase>(path, graph);
	}
	@Override
	public PathNode<DeBruijnNodeBase> createPathNode(PathNode<DeBruijnNodeBase> unsplit, int startIndex, int length, DeBruijnGraphBase<DeBruijnNodeBase> graph) {
		return new PathNode<DeBruijnNodeBase>(unsplit, startIndex, length, graph);
	}
	@Override
	public PathNode<DeBruijnNodeBase> createPathNode(DeBruijnGraphBase<DeBruijnNodeBase> graph, Iterable<PathNode<DeBruijnNodeBase>> path) {
		return new PathNode<DeBruijnNodeBase>(graph, path);
	}
}
