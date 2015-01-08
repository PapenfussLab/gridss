package au.edu.wehi.idsv;


public class SoftClipParameters {
	/**
	 * Minimum soft clip length to be considered evidence
	 */
	public int minLength = 4;
	/**
	 * Minimum MAPQ of read to considered evidence
	 */
	public int minReadMapq = 5;
	/**
	 * Minimum anchor percent identity to considered evidence
	 * 0-100
	 */
	public float minAnchorIdentity = 95;
	public AdapterHelper adapters = new AdapterHelper(new String[] {
			"AGATCGGAAGAG", // Illumina Universal Adapter
			"ATGGAATTCTCG", // Illumina Small RNA Adapter
			"CTGTCTCTTATA", // Nextera Transposase Sequence
	});
}
