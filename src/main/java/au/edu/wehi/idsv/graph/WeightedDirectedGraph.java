package au.edu.wehi.idsv.graph;


public interface WeightedDirectedGraph<T> extends DirectedGraph<T> {
	int getWeight(T node);
}
