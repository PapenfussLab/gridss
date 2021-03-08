package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class DownsamplingConfiguration {
	public static final String CONFIGURATION_PREFIX = "downsample";
	public DownsamplingConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		acceptDensityPortion = config.getDouble("acceptDensityPortion");
		targetEvidenceDensity = config.getDouble("targetEvidenceDensity");
		minimumDensityWindowSize = config.getInt("minimumDensityWindowSize");
		densityDownsampleRateClippedReads = config.getFloat("densityDownsampleRateClippedReads");
		densityDownsampleRateDiscordantReads = config.getFloat("densityDownsampleRateDiscordantReads");
	}
	/**
	 * Evidence per base to assemble without filtering
	 */
	public double acceptDensityPortion;
	/**
	 * Filter evidence to expect this maximum average evidence per base within window 
	 */
	public double targetEvidenceDensity;
	/**
	 * Minimum window size for density calculation 
	 */
	public int minimumDensityWindowSize;
	/**
	 * Downsampling rate for soft clipped, split read and indel evidence
	 */
	public final float densityDownsampleRateClippedReads;
	/**
	 * Downsampling rate for read pair evidence
	 */
	public final float densityDownsampleRateDiscordantReads;
}