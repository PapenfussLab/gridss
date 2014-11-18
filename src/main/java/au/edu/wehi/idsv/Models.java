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
	private static final Log log = Log.getInstance(AssemblyFactory.class);
	/**
	 * Calculates a somatic p-value for the given call
	 * Null hypothesis: germline and somatic BAFs are the same
	 * @param context processing context 
	 * @param call structural variant
	 * @return somatic p-value
	 */
	public static double somaticPvalue(ProcessingContext context, VariantContextDirectedEvidence call) {
		int variantNormal = 0, variantTumour = 0;
		int refNormal = call.getReferenceReadCount(EvidenceSubset.NORMAL) + call.getReferenceReadPairCount(EvidenceSubset.NORMAL);
		int refTumour = call.getReferenceReadCount(EvidenceSubset.TUMOUR) + call.getReferenceReadPairCount(EvidenceSubset.TUMOUR);
		if (call.getMappedEvidenceCountAssembly() > 0) {
			// just use the assembled evidence as we know this supports our variant 
			variantNormal = call.getAssemblySupportCount(EvidenceSubset.NORMAL);
			variantTumour = call.getAssemblySupportCount(EvidenceSubset.TUMOUR);
		} else {
			variantNormal = call.getEvidenceCount(EvidenceSubset.NORMAL);
			variantTumour = call.getEvidenceCount(EvidenceSubset.TUMOUR);
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
		if (total == 0) return 1.0d; // TODO: calculate p-value from both sides of the breakend
		double expectedBaf = variantTotal / total;
		BinomialDistribution tumourDist = new BinomialDistribution(tumourTotal, expectedBaf);
		// Single sided p-value: variant in tumour and not in normal. TODO: somatic LOH
		// TODO: is this numerically stable and performant for large n?
		// TODO: should we use a normal approximation for large N? 
		double somaticP = tumourDist.cumulativeProbability(variantTumour - 1, tumourTotal); // lower bound is exclusive
		return somaticP;
	}
	/**
	 * Log-likelihood ratio of existence of a structural variation supporting allele vs all reference alleles
	 * @return Log-likelihood ratio
	 */
	public static double llr(DirectedEvidence e) {
		if (e == null) throw new NullPointerException();
		if (e instanceof RealignedSoftClipEvidence) {
			return rscLlr((RealignedSoftClipEvidence)e);
		} else if (e instanceof SoftClipEvidence) {
			return scLlr((SoftClipEvidence)e);
		} else if (e instanceof DiscordantReadPair) {
			return dpLlr((DiscordantReadPair)e);
		} else if (e instanceof UnmappedMateReadPair) {
			return oeaLlr((UnmappedMateReadPair)e);
		} else if (e instanceof AssemblyEvidence) {
			return assemblyLlr((AssemblyEvidence)e);
		} else if (e instanceof VariantContextDirectedEvidence) {
			return ((VariantContextDirectedEvidence)e).getPhredScaledQual();
		}
		throw new IllegalArgumentException("Unknown evidence type " + e.getClass().toString());
	}
	private static double assemblyLlr(AssemblyEvidence e) {
		int evidenceCount = e.getAssemblySupportCountReadPair(EvidenceSubset.ALL) + e.getAssemblySupportCountSoftClip(EvidenceSubset.ALL);
		if (e instanceof DirectedBreakpoint) {
			return ((DirectedBreakpoint)e).getRemoteMapq() * evidenceCount;
		} else {
			return e.getLocalMapq() * evidenceCount;
		}
	}
	private static double oeaLlr(UnmappedMateReadPair e) {
		return Math.min(15, e.getLocalMapq());
	}
	private static double dpLlr(DiscordantReadPair e) {
		return Math.min(15, Math.min(e.getLocalMapq(), e.getRemoteMapq()));
	}
	private static double scLlr(SoftClipEvidence e) {
		return Math.min(Math.min(10, e.getSoftClipLength()), e.getLocalMapq());
	}
	private static double rscLlr(RealignedSoftClipEvidence e) {
		return Math.min(20, Math.min(e.getRemoteMapq(), e.getLocalMapq()));
	}
	// TODO: proper 95% Confidence Interval instead of hard limits on the bounds
	/**
	 * @param evidence
	 * @return
	 */
	public static BreakendSummary calculateBreakend(List<? extends DirectedEvidence> evidence) {
		if (evidence == null || evidence.size() == 0) throw new IllegalArgumentException("No evidence supplied");
		Collections.sort(evidence, DirectedEvidence.ByLlrDesc);
		// Start with best evidence
		BreakendSummary bounds = evidence.get(0).getBreakendSummary();
		bounds = new BreakendSummary(bounds.referenceIndex, bounds.direction, bounds.start, bounds.end);
		for (DirectedEvidence e : evidence) {
			// Reduce as we can
			BreakendSummary newbounds = BreakendSummary.overlapOf(bounds, e.getBreakendSummary());
			if (newbounds == null) {
				log.debug("Inconsisent support set");
			} else {
				bounds = newbounds;
			}
		}
		return bounds;
	}
}
