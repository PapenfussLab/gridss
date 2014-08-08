package au.edu.wehi.idsv;

public class ReadPairParameters {
	/**
	 * Use SAM flag to extract read pairs
	 */
	public boolean useProperPairFlag = true;
	/**
	 * Percentage (0.0-1.0) of reads considered concordant and ignored for read pair analysis.
	 */
	public float concordantPercent = 0;
}
