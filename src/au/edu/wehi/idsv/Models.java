package au.edu.wehi.idsv;

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
	 * @param call
	 * @return
	 */
	public static double somaticp(ProcessingContext context, VariantContextDirectedEvidence call) {
		int variantNormal = 0, variantTumour = 0;
		int refNormal = call.getReferenceReadCount(Subset.NORMAL) + call.getReferenceReadPairCount(Subset.NORMAL);
		int refTumour = call.getReferenceReadCount(Subset.TUMOUR) + call.getReferenceReadPairCount(Subset.TUMOUR);
		if (call.getMappedEvidenceCountAssembly() > 0) {
			// just use the assembly
			variantNormal = call.getAssemblySupportCount(Subset.NORMAL);
			variantTumour = call.getAssemblySupportCount(Subset.TUMOUR);
		} else {
			variantNormal = call.getEvidenceCount(Subset.NORMAL);
			variantTumour = call.getEvidenceCount(Subset.TUMOUR);
		}
		// TODO:
		// null = germline variant
		// hypothesis = somatic variant
		return variantTumour > 10 * variantNormal ? 1 : 0;
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
