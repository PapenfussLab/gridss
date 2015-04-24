package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.graph.WeightedDirectedGraph;

public interface DeBruijnGraph<T> extends WeightedDirectedGraph<T> {
	int getK();
	long getKmer(T node);
	boolean isReference(T node);
}
