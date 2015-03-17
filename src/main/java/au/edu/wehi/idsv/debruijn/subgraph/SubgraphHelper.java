package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

class SubgraphHelper {
	public static int getScore(
			List<SubgraphPathNode> current,
			Function<SubgraphPathNode, Integer> perNodeScoringFunction,
			boolean scoreReference,
			boolean scoreNonReference) {
		int sum = 0;
		for (SubgraphPathNode n : current) {
			sum += getScore(n, perNodeScoringFunction, scoreReference, scoreNonReference);
		}
		return sum;
	}
	public static int getScore(
			SubgraphPathNode node,
			Function<SubgraphPathNode, Integer> scoringFunction,
			boolean scoreReference,
			boolean scoreNonReference) {
		if ((node.containsNonReferenceKmer() && scoreNonReference) ||
			 (node.containsReferenceKmer() && scoreReference)) {
			return scoringFunction.apply(node);
		}
		return 0;
	}
	public static Ordering<SubgraphPathNode> getScoreOrderDesc(final Function<SubgraphPathNode, Integer> scoringFunction) {
		return new Ordering<SubgraphPathNode>() {
			@Override
			public int compare(SubgraphPathNode arg0, SubgraphPathNode arg1) {
				return Ints.compare(scoringFunction.apply(arg1), scoringFunction.apply(arg0));
			}};
	}
	public static final Function<SubgraphPathNode, Integer> SCORE_TOTAL_KMER = new Function<SubgraphPathNode, Integer>() {
		public Integer apply(SubgraphPathNode arg) {
			return arg.getWeight();
		}
	};
	public static final Function<SubgraphPathNode, Integer> SCORE_MAX_KMER = new Function<SubgraphPathNode, Integer>() {
		public Integer apply(SubgraphPathNode arg) {
			return arg.getMaxKmerWeight();
		}
	};
}
