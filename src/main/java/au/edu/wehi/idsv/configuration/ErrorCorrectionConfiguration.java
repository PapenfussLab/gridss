package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class ErrorCorrectionConfiguration {
	public static final String CONFIGURATION_PREFIX = "errorCorrection";
	public ErrorCorrectionConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		kmerErrorCorrectionMultiple = config.getFloat("kmerErrorCorrectionMultiple");
		k = config.getInt("k");
		maxCorrectionsInKmer = config.getInt("maxCorrectionsInKmer");
		deduplicateReadKmers = config.getBoolean("deduplicateReadKmers");
	}
	/**
	 * Extent to which an adjacent kmer should be more supported before error correction
	 */
	public float kmerErrorCorrectionMultiple;
	/**
	 * error correction kmer size
	 */
	public int k;
	/**
	 * Maximum number of corrections in a kmer before abandoning error correction
	 */
	public int maxCorrectionsInKmer;
	/**
	 * Count kmers once per read, or once per occurrence.
	 */
	public boolean deduplicateReadKmers;
}