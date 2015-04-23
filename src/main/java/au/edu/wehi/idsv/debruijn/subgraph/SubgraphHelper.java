package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

class SubgraphHelper {
	public static int getScore(
			List<PathNode> current,
			Function<PathNode, Integer> perNodeScoringFunction,
			boolean scoreReference,
			boolean scoreNonReference) {
		int sum = 0;
		for (PathNode n : current) {
			sum += getScore(n, perNodeScoringFunction, scoreReference, scoreNonReference);
		}
		return sum;
	}
	public static int getScore(
			PathNode node,
			Function<PathNode, Integer> scoringFunction,
			boolean scoreReference,
			boolean scoreNonReference) {
		if ((node.containsNonReferenceKmer() && scoreNonReference) ||
			 (node.containsReferenceKmer() && scoreReference)) {
			return scoringFunction.apply(node);
		}
		return 0;
	}
	public static Ordering<PathNode> getScoreOrderDesc(final Function<PathNode, Integer> scoringFunction) {
		return new Ordering<PathNode>() {
			@Override
			public int compare(PathNode arg0, PathNode arg1) {
				return Ints.compare(scoringFunction.apply(arg1), scoringFunction.apply(arg0));
			}};
	}
	public static final Function<PathNode, Integer> SCORE_TOTAL_KMER = new Function<PathNode, Integer>() {
		public Integer apply(PathNode arg) {
			return arg.weight();
		}
	};
	public static final Function<PathNode, Integer> SCORE_MAX_KMER = new Function<PathNode, Integer>() {
		public Integer apply(PathNode arg) {
			return arg.getMaxKmerWeight();
		}
	};
}
