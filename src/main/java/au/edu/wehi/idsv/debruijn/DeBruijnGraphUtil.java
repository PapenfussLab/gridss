package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.graph.WeightedDirectedGraph;

public class DeBruijnGraphUtil {
	public static <T> int score(WeightedDirectedGraph<T> graph, T n) {
		return graph.getWeight(n);
	}
	public static <T> int scoreAll(WeightedDirectedGraph<T> graph, Iterable<? extends T> path) {
		int sum = 0;
		for (T n : path) {
			sum += score(graph, n);
		}
		return sum;
	}
	public static <PN extends PathNode<?>> int score(DeBruijnGraph<PN> graph, PN node, boolean scoreReference, boolean scoreNonReference) {
		boolean isRef = graph.isReference(node);
		if ((isRef && scoreReference) || (!isRef && scoreNonReference)) {
			// scoring function is the node weight
			return node.weight();
		}
		return 0;
	}
	public static <PN extends PathNode<?>> int scorePath(DeBruijnGraph<PN> graph, Iterable<? extends PN> path, boolean scoreReference, boolean scoreNonReference) {
		int sum = 0;
		for (PN n : path) {
			sum += score(graph, n, scoreReference, scoreNonReference);
		}
		return sum;
	}
}
