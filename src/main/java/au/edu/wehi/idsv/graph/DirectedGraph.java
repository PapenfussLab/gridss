package au.edu.wehi.idsv.graph;

import java.util.Collection;
import java.util.List;

public interface DirectedGraph<T> {
	List<T> next(T node);
	List<T> prev(T node);
	void removeNode(T node);
	void removeEdge(T source, T sink);
	void addNode(T node);
	void addEdge(T source, T sink);
	Collection<T> allNodes();
	String toString(Iterable<? extends T> path);
}
