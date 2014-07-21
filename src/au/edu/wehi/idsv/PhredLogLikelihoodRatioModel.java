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
			if (e.isSimpleAssembly()) {
				return assembly((VariantContextDirectedEvidence)e);
			}
			return llr((VariantContextDirectedEvidence)e);
		}
		throw new IllegalArgumentException("Unknown evidence type " + e.class.toString());
	}
	private static void assembly(VariantContextDirectedEvidence e) {
		// TODO Auto-generated method stub
		return 0;
	}
	private static float oea(UnmappedMateReadPair e) {
		Math.min(15, e.getLocalMapq());
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
