package au.edu.wehi.idsv;

public class SoftClipParameters {
	/**
	 * Minimum soft clip length to be considered evidence
	 */
	public int minLength = 4;
	/**
	 * Minimum MAPQ for realignment
	 */
	public int maxMapq = 5;
	/**
	 * Minimum anchor percent identity
	 * 0-100
	 */
	public float minAnchorIdentity = 95;
}
