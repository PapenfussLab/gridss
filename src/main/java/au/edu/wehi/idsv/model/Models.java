package au.edu.wehi.idsv.model;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.graph.MaximumCliqueIntervalGraph;
import au.edu.wehi.idsv.graph.MaximumCliqueIntervalGraph.Node;
import au.edu.wehi.idsv.graph.ScalingHelper;
import au.edu.wehi.idsv.util.MathUtil;

/**
 * variant/reference Log-likelihood statistical model
 * 
 * @author Daniel Cameron
 *
 */
public class Models {
	//private static final Log log = Log.getInstance(Models.class);
	/**
	 * Calculates the most likely breakend interval for the given evidence 
	 * @param evidence
	 * @return breakend interval with highest total evidence quality
	 */
	public static BreakendSummary calculateBreakend(LinearGenomicCoordinate lgc, List<? extends DirectedEvidence> evidence) {
		if (evidence == null || evidence.size() == 0) throw new IllegalArgumentException("No evidence supplied");
		return calculateBreakend(lgc, Lists.transform(evidence, new Function<DirectedEvidence, BreakendSummary>() {
				@Override
				public BreakendSummary apply(DirectedEvidence input) {
					return input.getBreakendSummary();
				}
			}), Lists.transform(evidence, new Function<DirectedEvidence, Long>() {
				@Override
				public Long apply(DirectedEvidence input) {
					return ScalingHelper.toScaledWeight(input.getBreakendQual());
				}
			}));
	}
	/**
	 * Calculates the most likely breakend interval for the given evidence 
	 * @param evidence
	 * @return breakend interval with highest total evidence quality
	 */
	public static BreakendSummary calculateBreakend(LinearGenomicCoordinate lgc, List<BreakendSummary> bs, List<Long> weights) {
		if (bs == null || bs.size() == 0) throw new IllegalArgumentException("No evidence supplied");
		if (weights.size() != bs.size()) throw new IllegalArgumentException("Array lenght mismatch");
		Node fwd = maximalInterval(lgc, BreakendDirection.Forward, bs, weights);
		Node bwd = maximalInterval(lgc, BreakendDirection.Backward, bs, weights);
		if (fwd == null && bwd == null) {
			// all evidence is insignificant, just return something as we're going to get filtered anyway
			BreakendSummary fallback = bs.get(0);
			if (fallback instanceof BreakpointSummary) {
				BreakpointSummary bp = (BreakpointSummary)fallback;
				fallback = bp.localBreakend();
			}
			return fallback;
		}
		Node node = fwd;
		BreakendDirection dir = BreakendDirection.Forward;
		if (fwd == null || (bwd != null && fwd.weight < bwd.weight)) {
			node = bwd;
			dir = BreakendDirection.Backward;
		}
		assert(lgc.getReferenceIndex(node.start) == lgc.getReferenceIndex(node.stop));
		int start = lgc.getReferencePosition(node.start);
		int end = lgc.getReferencePosition(node.stop);
		return new BreakendSummary(lgc.getReferenceIndex(node.start), dir, MathUtil.average(start, end), start, end);
	}
	private static Node maximalInterval(LinearGenomicCoordinate lgc, BreakendDirection dir, List<BreakendSummary> breaks, List<Long> weights) {
		MaximumCliqueIntervalGraph calc = new MaximumCliqueIntervalGraph();
		List<Node> nodes = new ArrayList<Node>(breaks.size());
		for (int i = 0; i < breaks.size(); i++) {
			BreakendSummary bs = breaks.get(i);
			long weight = weights.get(i);
			if (bs == null) continue;
			if (bs.direction != dir) continue;
			if (weight > 0) {
				nodes.add(new Node(
						lgc.getLinearCoordinate(bs.referenceIndex, bs.start), 
						lgc.getLinearCoordinate(bs.referenceIndex, bs.end),
						weight));
			}
		}
		if (nodes.size() == 0) return null;
		return calc.calculateMaximumClique(nodes);
	}
}
