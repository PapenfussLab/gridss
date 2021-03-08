package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class ErrorCorrectionConfiguration {
	public static final String CONFIGURATION_PREFIX = "errorCorrection";
	public ErrorCorrectionConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		maxBaseMismatchForCollapse = config.getInt("maxBaseMismatchForCollapse");
		collapseBubblesOnly = config.getBoolean("collapseBubblesOnly");
		maxPathCollapseLengthMultiple = config.getFloat("maxPathCollapseLengthMultiple");
		kmerErrorCorrectionMultiple = config.getFloat("kmerErrorCorrectionMultiple");
		k = config.getInt("k");
	}
	/**
	 * Maximum of base mismatches for de bruijn kmer paths to be merged   
	 */
	public int maxBaseMismatchForCollapse;
	/**
	 * Only collapse bubble path with a single entry and exit kmer choice
	 */
	public boolean collapseBubblesOnly;
	/**
	 * Maximum length of path to collapse
	 * Units are multiples of max support width (ie largest max fragment size)
	 */
	public float maxPathCollapseLengthMultiple;
	/**
	 * Extent to which an adjacent kmer should be more supported before error correction
	 */
	public float kmerErrorCorrectionMultiple;
	/**
	 * error correction kmer size
	 */
	public int k;
	public int maxPathCollapseLengthInBases(int readLength) { return (int)(maxPathCollapseLengthMultiple * readLength); }
}