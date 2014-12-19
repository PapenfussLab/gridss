package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;

import au.edu.wehi.idsv.metrics.IdsvMetrics;
import au.edu.wehi.idsv.metrics.InsertSizeDistribution;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

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
		double llr = -1;
		if (e instanceof RealignedRemoteSoftClipEvidence) {
			llr = rrscLlr((RealignedRemoteSoftClipEvidence)e);
		} else if (e instanceof RealignedSoftClipEvidence) {
			llr = rscLlr((RealignedSoftClipEvidence)e);
		} else if (e instanceof SoftClipEvidence) {
			llr = scLlr((SoftClipEvidence)e);
		} else if (e instanceof DiscordantReadPair) {
			llr = dpLlr((DiscordantReadPair)e);
		} else if (e instanceof UnmappedMateReadPair) {
			llr = oeaLlr((UnmappedMateReadPair)e);
		} else if (e instanceof AssemblyEvidence) {
			llr = assemblyLlr((AssemblyEvidence)e);
		} else if (e instanceof VariantContextDirectedEvidence) {
			llr = ((VariantContextDirectedEvidence)e).getPhredScaledQual();
		} else {
			throw new IllegalArgumentException("Unknown evidence type " + e.getClass().toString());
		}
		if (llr < 0) throw new RuntimeException("Sanity check failure: negative llr");
		return llr;
		
	}
	private static double assemblyLlr(AssemblyEvidence e) {
		int evidenceCount = e.getAssemblySupportCountReadPair(EvidenceSubset.ALL) + e.getAssemblySupportCountSoftClip(EvidenceSubset.ALL);
		if (e instanceof DirectedBreakpoint) {
			return Math.max(1, ((DirectedBreakpoint)e).getRemoteMapq()) * evidenceCount;
		} else {
			return Math.max(1, e.getLocalMapq()) * evidenceCount;
		}
	}
	private static double oeaLlr(UnmappedMateReadPair e) {
		return oeaPhred(e.getEvidenceSource(), e.getLocalMapq());
	}
	private static double scLlr(SoftClipEvidence e) {
		return scPhred(e.getEvidenceSource(), e.getSoftClipLength(), e.getLocalMapq());
	}
	private static double rscLlr(RealignedSoftClipEvidence e) {
		return scPhred(e.getEvidenceSource(), e.getSoftClipLength(), e.getLocalMapq(), e.getRemoteMapq());
	}
	private static double rrscLlr(RealignedRemoteSoftClipEvidence e) {
		return scPhred(e.getEvidenceSource(), e.getOriginalSoftClipLength(), e.getRemoteMapq(), e.getLocalMapq());
	}
	private static double scPhred(SAMEvidenceSource source, int clipLength, int localMapq, int remoteMapq) {
		double score = scPhred(source, clipLength, localMapq);
		score = Math.min(score, remoteMapq);
		return score;
	}
	public static int getAssemblyScore(Collection<DirectedEvidence> support) {
		int score = 0;
		for (DirectedEvidence e : support) {
			score += llr(e);
		}
		return score;
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
		Collections.sort(evidence, DirectedEvidence.ByLlrDesc);
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
