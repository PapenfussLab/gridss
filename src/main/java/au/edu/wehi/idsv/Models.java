package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;

import au.edu.wehi.idsv.graph.MaximumCliqueIntervalGraph;
import au.edu.wehi.idsv.graph.MaximumCliqueIntervalGraph.Node;
import au.edu.wehi.idsv.graph.ScalingHelper;

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
		if (evidence == null || evidence.size() == 0) {
			throw new IllegalArgumentException("No evidence supplied");
		}
		BreakendDirection direction = evidence.get(0).getBreakendSummary().direction;
		MaximumCliqueIntervalGraph calc = new MaximumCliqueIntervalGraph();
		List<Node> nodes = new ArrayList<Node>(evidence.size());
		for (DirectedEvidence e : evidence) {
			long weight = ScalingHelper.toScaledWeight(e.getBreakendQual());
			if (weight > 0) {
			BreakendSummary bs = e.getBreakendSummary();
			assert(bs.direction == direction);
			nodes.add(new Node(
					lgc.getLinearCoordinate(bs.referenceIndex, bs.start), 
					lgc.getLinearCoordinate(bs.referenceIndex, bs.end),
					weight));
			}
		}
		Node call = calc.calculateMaximumClique(nodes);
		return new BreakendSummary(lgc.getReferenceIndex(call.start), direction, lgc.getReferencePosition(call.start), lgc.getReferencePosition(call.stop));
	}
}
