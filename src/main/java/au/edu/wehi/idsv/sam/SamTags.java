package au.edu.wehi.idsv.sam;

import java.util.List;

import com.google.common.collect.ImmutableList;

public class SamTags {
	/**
	 * Multiple non-chimeric mapping locations have been reported by the aligner
	 * for at least one read in his fragments.  This fields contains the number
	 * of additional mappings for any read in the fragment.
	 * 
	 * For example, if the reads in a read pair have 5 and 10 mapping locations
	 * respectively, mm will be set to 4+9=13.
	 */
	public static final String MULTIMAPPING_FRAGMENT = "mm";
	/**
	 * Contig of the remainder of the realigned read
	 */
	public static final String REALIGNMENT_REFERENCE_INDEX = "rr";
	/**
	 * Position of the remainder of the realigned read
	 */
	public static final String REALIGNMENT_POSITION = "rp";
	/**
	 * Filters applied to assembly
	 */
	public static final String ASSEMBLY_FILTERS = "af";
	/**
	 * EvidenceID of assembly components
	 */
	public static final String EVIDENCEID = "ez";
	/**
	 * Read names of supporting fragments.
	 * 
	 * Per-category encoding is achieved via consecutive separators
	 */
	public static final String ASSEMBLY_SUPPORTING_FRAGMENTS = "sf";
	public static final String ASSEMBLY_DIRECTION = "ad";
	public static final String UNANCHORED = "ua";
	// Per category aggregations
	public static final String ASSEMBLY_READPAIR_COUNT = "dc";
	public static final String ASSEMBLY_SOFTCLIP_COUNT = "sc";
	public static final String ASSEMBLY_STRAND_BIAS = "sb";
	public static final String ASSEMBLY_READPAIR_QUAL = "dq";
	public static final String ASSEMBLY_SOFTCLIP_QUAL = "sq";
	public static final String ASSEMBLY_READPAIR_LENGTH_MAX = "dl";
	/***
	 * Pair of integers indicating how many bases at the start/end of the contig were truncated
	 * to ensure the contig anchor did not overrun the contig bounds
	 */
	public static final String ASSEMBLY_ANCHOR_TRUNCATION = "at";
	public static final List<String> ASSEMBLY_ANNOTATIONS = ImmutableList.of(
			ASSEMBLY_FILTERS,
			EVIDENCEID,
			ASSEMBLY_DIRECTION,
			UNANCHORED,
			ASSEMBLY_READPAIR_COUNT,
			ASSEMBLY_SOFTCLIP_COUNT,
			ASSEMBLY_READPAIR_QUAL,
			ASSEMBLY_SOFTCLIP_QUAL,
			ASSEMBLY_READPAIR_LENGTH_MAX,
			ASSEMBLY_STRAND_BIAS,
			ASSEMBLY_ANCHOR_TRUNCATION
			);
}
