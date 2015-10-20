package au.edu.wehi.idsv.configuration;

import au.edu.wehi.idsv.AdapterHelper;


public class SoftClipConfiguration {
	/**
	 * Minimum average breakend quality score to be considered a valid soft clip
	 * This filters out reads that are soft-clipped due to sequencing errors  
	 */
	public float minAverageQual = 5;
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
	/**
	 * Minimum entropy of anchored sequence (in bits) (Shannon entropy)
	 */
	public double minAnchorEntropy = 0.5;
	public AdapterHelper adapters = new AdapterHelper(new String[] {
			"AGATCGGAAGAG", // Illumina Universal Adapter
			"ATGGAATTCTCG", // Illumina Small RNA Adapter
			"CTGTCTCTTATA", // Nextera Transposase Sequence
	});
}
