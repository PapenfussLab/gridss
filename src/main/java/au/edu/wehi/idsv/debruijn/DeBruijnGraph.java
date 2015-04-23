package au.edu.wehi.idsv.debruijn;

import java.util.Iterator;

import au.edu.wehi.idsv.graph.WeightedDirectedGraph;

public interface DeBruijnGraph<T> extends WeightedDirectedGraph<T> {
	long getKmer(T node);
	boolean isReference(T node);
	int basesDifferent(Iterator<T> pathA, Iterator<T> pathB);
}
