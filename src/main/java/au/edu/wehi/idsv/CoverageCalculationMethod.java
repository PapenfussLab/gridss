package au.edu.wehi.idsv;

public enum CoverageCalculationMethod {
	/**
	 * Calculate physical coverage by using the inferring fragment size assuming no structural variation
	 */
	FRAGMENT,
	/**
	 * Calculate read coverage
	 */
	READ,
}
