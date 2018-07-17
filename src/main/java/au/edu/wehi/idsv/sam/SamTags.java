package au.edu.wehi.idsv.sam;

public class SamTags {
	public static final String ASSEMBLY_DIRECTION = "ad";
	public static final String UNANCHORED = "ua";
	public static final String ASSEMBLY_STRAND_BIAS = "sb";
	public static final String ASSEMBLY_EVIDENCE_OFFSET_START = "os";
	public static final String ASSEMBLY_EVIDENCE_OFFSET_END = "oe";
	public static final String ASSEMBLY_EVIDENCE_TYPE = "et";
	public static final String ASSEMBLY_EVIDENCE_CATEGORY = "ec";
	public static final String ASSEMBLY_EVIDENCE_QUAL = "eq";
	public static final String ASSEMBLY_EVIDENCE_EVIDENCEID = "ez";
	public static final String ASSEMBLY_EVIDENCE_FRAGMENTID = "ef";
	public static final String ASSEMBLY_MAX_READ_LENGTH = "dl";
	/***
	 * Pair of integers indicating how many bases at the start/end of the contig were truncated
	 * to ensure the contig anchor did not overrun the contig bounds
	 */
	public static final String ASSEMBLY_ANCHOR_TRUNCATION = "at";
}
