package au.edu.wehi.idsv;

public class SoftClipParameters {
	/**
	 * Minimum soft clip length to be considered evidence
	 */
	public int minLength = 4;
	/**
	 * Minimum MAPQ of read to consider for realignment
	 */
	public int minReadMapq = 5;
	/**
	 * Minimum anchor percent identity to consider for realignment
	 * 0-100
	 */
	public float minAnchorIdentity = 95;
	public boolean meetsEvidenceCritera(SoftClipEvidence sce) {
		return sce.getMappingQuality() >= minReadMapq
				&& sce.getSoftClipLength() >= minLength
				&& sce.getAlignedPercentIdentity() >= minAnchorIdentity;
	}
}
