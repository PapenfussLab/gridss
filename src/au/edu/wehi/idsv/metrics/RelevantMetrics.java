package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.SamPairUtil.PairOrientation;

public interface RelevantMetrics {

	/**
	 * Gets the median fragment size
	 * @return median fragment size
	 */
	public abstract double getMedianFragmentSize();

	/**
	 * Gets the standard deviation of the fragment size
	 * @return fragment size standard deviation
	 */
	public abstract double getFragmentSizeStdDev();

	/**
	 * Gets the maximum expected fragment size
	 * @return longest expected fragment size
	 */
	public abstract int getMaxFragmentSize();
	/**
	 * Maximum read length
	 * @return longest read length
	 */
	public abstract int getMaxReadLength();
	public abstract PairOrientation getPairOrientation();

}