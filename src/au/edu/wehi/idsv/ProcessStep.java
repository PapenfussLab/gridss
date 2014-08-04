package au.edu.wehi.idsv;

import java.util.EnumSet;

/**
 * Data transformation steps performed by IDSV  
 * 
 * @author Daniel Cameron
 *
 */
public enum ProcessStep {
	/**
	 * Calculates relevant metrics from the input SAM/BAM
	 */
	CALCULATE_METRICS,
	/**
	 * Extracts soft clip evidence from the input SAM/BAM
	 */
	EXTRACT_SOFT_CLIPS,
	/**
	 * Extracts non-concordant read pairs from the input SAM/BAM
	 */
	EXTRACT_READ_PAIRS,
	/**
	 * Extracts read pair mates from the input SAM/BAM
	 */
	EXTRACT_READ_MATES,
	/**
	 * Sorts read pair mates according to partner position 
	 */
	SORT_READ_MATES,
	/**
	 * Generates breakends assemblies from the soft clip and
	 * read pair evidence
	 */
	ASSEMBLE_BREAKENDS,
	/**
	 * Aligns soft clipped sequences to the reference
	 */
	REALIGN_SOFT_CLIPS,
	/**
	 * Aligns breakend assemblies to the reference
	 */
	REALIGN_ASSEMLIES,
	/**
	 * Sorts soft clips by realigned position 
	 */
	SORT_REALIGNED_SOFT_CLIPS,
	/**
	 * Sorts assembles by realigned breakend position 
	 */
	SORT_REALIGNED_ASSEMBLIES,
	/**
	 * Calls putative structural variants from the
	 * avaialble evidence 
	 */
	CALL_STRUCTURAL_VARIANTS,
	/**
	 * Annotates variants
	 */
	ANNOTATE_VARIANTS;
	/**
	 * All processing steps
	 */
	public static final EnumSet<ProcessStep> ALL_STEPS = EnumSet.allOf(ProcessStep.class);
}
