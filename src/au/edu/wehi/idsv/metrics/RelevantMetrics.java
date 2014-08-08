package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.SamPairUtil.PairOrientation;

public interface RelevantMetrics {

	/**
	 * Gets the median fragment size
	 * @return median fragment size
	 */
	double getMedianFragmentSize();

	/**
	 * Gets the standard deviation of the fragment size
	 * @return fragment size standard deviation
	 */
	double getFragmentSizeStdDev();

	/**
	 * Gets the maximum expected fragment size
	 * @return longest expected fragment size
	 */
	int getMaxFragmentSize();
	/**
	 * Maximum read length
	 * @return longest read length
	 */
	int getMaxReadLength();
	PairOrientation getPairOrientation();
	/**
	 * Gets the distribution of insert sizes
	 * @return insert size distribution
	 */
	InsertSizeDistribution getInsertSizeDistribution();
}