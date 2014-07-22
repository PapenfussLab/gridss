package au.edu.wehi.idsv;

/**
 * variant/reference Log-likelihood statistical model
 * 
 * @author cameron.d
 *
 */
public class PhredLogLikelihoodRatioModel {
	/**
	 * Log-likelihood ratio of existence of a structural variation supporting allele vs all reference alleles
	 * @return Log-likelihood ratio
	 */
	public static float llr(DirectedEvidence e) {
		if (e == null) throw new NullPointerException();
		if (e instanceof RealignedSoftClipEvidence) {
			return rsc((RealignedSoftClipEvidence)e);
		} else if (e instanceof SoftClipEvidence) {
			return sc((SoftClipEvidence)e);
		} else if (e instanceof DiscordantReadPair) {
			return dp((DiscordantReadPair)e);
		} else if (e instanceof UnmappedMateReadPair) {
			return oea((UnmappedMateReadPair)e);
		} else if (e instanceof VariantContextDirectedEvidence) {
			VariantContextDirectedEvidence vcde = (VariantContextDirectedEvidence)e;
			if (vcde.isSimpleAssembly()) {
				return assembly(vcde);
			}
			return llr((VariantContextDirectedEvidence)e);
		}
		throw new IllegalArgumentException("Unknown evidence type " + e.getClass().toString());
	}
	private static float assembly(VariantContextDirectedEvidence e) {
		if (e.getMappedEvidenceCountAssembly() > 0) {
			return Math.min(25, e.getMapqAssemblyRemoteMax()) * e.getAssemblySupportCount(null);
		} else {
			return 15 * e.getAssemblySupportCountReadPair(null) + 20 * e.getAssemblySupportCountSoftClip(null);
		}
	}
	private static float oea(UnmappedMateReadPair e) {
		return Math.min(15, e.getLocalMapq());
	}
	private static float dp(DiscordantReadPair e) {
		return Math.min(15, Math.min(e.getLocalMapq(), e.getRemoteMapq()));
	}
	private static float sc(SoftClipEvidence e) {
		return Math.min(Math.min(10, e.getSoftClipLength()), e.getLocalMapq());
	}
	private static float rsc(RealignedSoftClipEvidence e) {
		return Math.min(20, Math.min(e.getRemoteMapq(), e.getLocalMapq()));
	}
}
