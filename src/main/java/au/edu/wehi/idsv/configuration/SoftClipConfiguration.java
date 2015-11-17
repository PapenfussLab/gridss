package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class SoftClipConfiguration {
	public static final String CONFIGURATION_PREFIX = "softclip";
	public SoftClipConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		minAverageQual = config.getFloat("minAverageQual");
		minLength = config.getInt("minLength");
		minAnchorIdentity = config.getFloat("minAnchorIdentity");
	}
	/**
	 * Minimum average breakend quality score to be considered a valid soft clip
	 * This filters out reads that are soft-clipped due to sequencing errors  
	 */
	public float minAverageQual;
	/**
	 * Minimum soft clip length to be considered evidence
	 */
	public int minLength;
	/**
	 * Minimum anchor identity to considered evidence
	 */
	public float minAnchorIdentity;
}
