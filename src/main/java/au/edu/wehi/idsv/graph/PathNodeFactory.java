package au.edu.wehi.idsv.graph;

import java.util.Collection;

public interface PathNodeFactory<T, PN extends PathNode<T>> {
	PN createPathNode(Collection<T> path);
	PN splitPathNode(PN unsplit, int startIndex, int length);
	PN concatPathNodes(Iterable<PN> path);
}
