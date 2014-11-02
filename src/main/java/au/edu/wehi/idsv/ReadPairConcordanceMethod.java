package au.edu.wehi.idsv;

public enum ReadPairConcordanceMethod {
	/**
	 * Base read pair concordence based on proper pair flag
	 * assigned by aligner.
	 */
	SAM_FLAG,
	/**
	 * Base read pair concordance by taking the upper/lower bounds of read pairs
	 */
	PERCENTAGE,
	/**
	 * Base read pair concordance on fixed min/max bounds values 
	 */
	FIXED,
}
