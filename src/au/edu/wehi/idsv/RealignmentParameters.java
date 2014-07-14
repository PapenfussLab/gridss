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
	public boolean meetsCritera(SoftClipEvidence sce) {
		return sce.getSoftClipLength() >= minLength
				&& sce.getAverageClipQuality() >= minAverageQual;
	}
}
