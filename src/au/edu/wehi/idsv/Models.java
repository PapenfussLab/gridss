package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import org.apache.commons.math3.distribution.BinomialDistribution;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.vcf.VcfAttributes.Subset;

/**
 * variant/reference Log-likelihood statistical model
 * 
 * @author cameron.d
 *
 */
public class Models {
	/**
	 * Calculates a somatic p-value for the given call
	 * Null hypothesis: germline and somatic BAFs are the same
	 * @param context processing context 
	 * @param call structural variant
	 * @return somatic p-value
	 */
	public static double somaticPvalue(ProcessingContext context, VariantContextDirectedEvidence call) {
		int variantNormal = 0, variantTumour = 0;
		int refNormal = call.getReferenceReadCount(Subset.NORMAL) + call.getReferenceReadPairCount(Subset.NORMAL);
		int refTumour = call.getReferenceReadCount(Subset.TUMOUR) + call.getReferenceReadPairCount(Subset.TUMOUR);
		if (call.getMappedEvidenceCountAssembly() > 0) {
			// just use the assembled evidence as we know this supports our variant 
			variantNormal = call.getAssemblySupportCount(Subset.NORMAL);
			variantTumour = call.getAssemblySupportCount(Subset.TUMOUR);
		} else {
			variantNormal = call.getEvidenceCount(Subset.NORMAL);
			variantTumour = call.getEvidenceCount(Subset.TUMOUR);
		}
		assert(variantNormal >= 0);
		assert(variantTumour >= 0);
		assert(refNormal >= 0);
		assert(refTumour >= 0);
		int variantTotal = variantNormal + variantTumour;
		int refTotal = refNormal + refTumour;
		int normalTotal = variantNormal + refNormal;
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
		} else if (e instanceof VariantContextDirectedEvidence) {
			VariantContextDirectedEvidence vcde = (VariantContextDirectedEvidence)e;
			if (vcde.isSimpleAssembly()) {
				return assemblyLlr(vcde);
			}
			return (double)vcde.getPhredScaledQual();
		}
		throw new IllegalArgumentException("Unknown evidence type " + e.getClass().toString());
	}
	private static double assemblyLlr(VariantContextDirectedEvidence e) {
		if (e.getMappedEvidenceCountAssembly() > 0) {
			return Math.min(25, e.getMapqAssemblyRemoteMax()) * e.getAssemblySupportCount(null);
		} else {
			return 15 * e.getAssemblySupportCountReadPair(null) + 20 * e.getAssemblySupportCountSoftClip(null);
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
}
