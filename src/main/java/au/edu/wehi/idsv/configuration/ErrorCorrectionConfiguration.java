package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class ErrorCorrectionConfiguration {
	public static final String CONFIGURATION_PREFIX = "errorCorrection";
	public ErrorCorrectionConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		kmerErrorCorrectionMultiple = config.getFloat("kmerErrorCorrectionMultiple");
		k = config.getInt("k");
	}
	/**
	 * Extent to which an adjacent kmer should be more supported before error correction
	 */
	public float kmerErrorCorrectionMultiple;
	/**
	 * error correction kmer size
	 */
	public int k;
}