package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;

import au.edu.wehi.idsv.graph.MaximumCliqueIntervalGraph;
import au.edu.wehi.idsv.graph.MaximumCliqueIntervalGraph.Node;
import au.edu.wehi.idsv.graph.ScalingHelper;

import com.google.common.base.Function;
import com.google.common.collect.Lists;

/**
 * variant/reference Log-likelihood statistical model
 * 
 * @author cameron.d
 *
 */
public class Models {
	//private static final Log log = Log.getInstance(Models.class);
	/**
	 * Calculates a somatic p-value for the given call
	 * Null hypothesis: germline and somatic BAFs are the same
	 * @param context processing context 
	 * @param call structural variant
	 * @return somatic p-value
	 */
	public static float somaticPvalue(ProcessingContext context, VariantContextDirectedEvidence call) {
		int refNormal = call.getReferenceReadCount(EvidenceSubset.NORMAL) + call.getReferenceReadPairCount(EvidenceSubset.NORMAL);
		int refTumour = call.getReferenceReadCount(EvidenceSubset.TUMOUR) + call.getReferenceReadPairCount(EvidenceSubset.TUMOUR);
		int variantNormal = call.getBreakendEvidenceCountReadPair(EvidenceSubset.NORMAL) + call.getBreakendEvidenceCountSoftClip(EvidenceSubset.NORMAL);
		int variantTumour = call.getBreakendEvidenceCountReadPair(EvidenceSubset.TUMOUR) + call.getBreakendEvidenceCountSoftClip(EvidenceSubset.TUMOUR);
		if (call instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)call;
			variantNormal += bp.getBreakpointEvidenceCountReadPair(EvidenceSubset.NORMAL) + bp.getBreakpointEvidenceCountSoftClip(EvidenceSubset.NORMAL);
			variantTumour += bp.getBreakpointEvidenceCountReadPair(EvidenceSubset.TUMOUR) + bp.getBreakpointEvidenceCountSoftClip(EvidenceSubset.TUMOUR);
		}
		assert(variantNormal >= 0);
		assert(variantTumour >= 0);
		assert(refNormal >= 0);
		assert(refTumour >= 0);
		int variantTotal = variantNormal + variantTumour;
		//int refTotal = refNormal + refTumour;
		//int normalTotal = variantNormal + refNormal;
		int tumourTotal = variantTumour + refTumour;
		int total = variantNormal + variantTumour + refNormal + refTumour;
		if (total == 0) return 1; // TODO: calculate p-value from both sides of the breakend
		float expectedBaf = variantTotal / total;
		BinomialDistribution tumourDist = new BinomialDistribution(tumourTotal, expectedBaf);
		// Single sided p-value: variant in tumour and not in normal. TODO: somatic LOH
		// TODO: is this numerically stable and performant for large n?
		// TODO: should we use a normal approximation for large N? 
		float somaticP = (float)tumourDist.cumulativeProbability(variantTumour - 1, tumourTotal); // lower bound is exclusive
		return somaticP;
	}
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
		return new BreakendSummary(lgc.getReferenceIndex(node.start), dir, lgc.getReferencePosition(node.start), lgc.getReferencePosition(node.stop));
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
