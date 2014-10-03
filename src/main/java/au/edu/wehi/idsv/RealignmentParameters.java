package au.edu.wehi.idsv;


public class RealignmentParameters {
	/**
	 * Minimum breakend length to be considered for realignment
	 */
	public int minLength = 25;
	/**
	 * Minimum average breakend quality score to be considered for realignment 
	 */
	public float minAverageQual = 5;
	public boolean shouldRealignBreakend(SoftClipEvidence evidence) {
		if (evidence.getBreakendSummary() instanceof BreakpointSummary) return false;
		return evidence.getSoftClipLength() >= minLength
				&& evidence.getAverageClipQuality() >= minAverageQual;
	}
	public boolean shouldRealignBreakend(VariantContextDirectedEvidence evidence) {
		if (evidence.getBreakendSummary() instanceof BreakpointSummary) return false;
		return evidence.getBreakendSequence().length >= minLength;
				//&& evidence.getAssemblyQuality() >= minAverageQual;
	}
}
