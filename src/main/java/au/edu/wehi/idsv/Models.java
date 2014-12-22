package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;

/**
 * variant/reference Log-likelihood statistical model
 * 
 * @author cameron.d
 *
 */
public class Models {
	private static final Log log = Log.getInstance(Models.class);
	/**
	 * Calculates a somatic p-value for the given call
	 * Null hypothesis: germline and somatic BAFs are the same
	 * @param context processing context 
	 * @param call structural variant
	 * @return somatic p-value
	 */
	public static float somaticPvalue(ProcessingContext context, VariantContextDirectedEvidence call) {
		int variantNormal = 0, variantTumour = 0;
		int refNormal = call.getReferenceReadCount(EvidenceSubset.NORMAL) + call.getReferenceReadPairCount(EvidenceSubset.NORMAL);
		int refTumour = call.getReferenceReadCount(EvidenceSubset.TUMOUR) + call.getReferenceReadPairCount(EvidenceSubset.TUMOUR);
		if (call instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)call;
			variantNormal = bp.getBreakpointEvidenceCountReadPair(EvidenceSubset.NORMAL) + bp.getBreakpointEvidenceCountSoftClip(EvidenceSubset.NORMAL);
			variantTumour = bp.getBreakpointEvidenceCountReadPair(EvidenceSubset.NORMAL) + bp.getBreakpointEvidenceCountSoftClip(EvidenceSubset.NORMAL);

		} else {
			throw new RuntimeException("NYI");
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
	// TODO: proper 95% Confidence Interval instead of hard limits on the bounds
	/**
	 * @param evidence
	 * @return
	 */
	public static BreakendSummary calculateBreakend(List<? extends DirectedEvidence> evidence) {
		if (evidence == null || evidence.size() == 0) {
			throw new IllegalArgumentException("No evidence supplied");
		}
		boolean errorMessagePrinted = false;
		Collections.sort(evidence, DirectedEvidence.ByBreakendQualDesc);
		// Start with best evidence
		BreakendSummary bounds = evidence.get(0).getBreakendSummary();
		bounds = new BreakendSummary(bounds.referenceIndex, bounds.direction, bounds.start, bounds.end);
		for (DirectedEvidence e : evidence) {
			// Reduce as we can
			BreakendSummary newbounds = BreakendSummary.overlapOf(bounds, e.getBreakendSummary());
			if (newbounds == null) {
				if (!errorMessagePrinted) {
					printInconsisentSupportSetErrorMessage(evidence);
					errorMessagePrinted = true;
				}
			} else {
				bounds = newbounds;
			}
		}
		return bounds;
	}
	private static void printInconsisentSupportSetErrorMessage(List<? extends DirectedEvidence> evidence) {
		StringBuilder sb = new StringBuilder("Inconsisent support set {");
		for (DirectedEvidence e : evidence) {
			sb.append(e.getBreakendSummary().toString());
			sb.append(", ");
		}
		sb.append("}");
		log.debug(sb.toString());
	}
}
