package au.edu.wehi.idsv.sam;

import htsjdk.samtools.SAMTag;

public class SamTags {
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
	public static final String ASSEMLBY_FILTERS = "af";
	/**
	 * EvidenceID of assembly components
	 */
	public static final String EVIDENCEID = "ez";
	public static final String ASSEMBLY_DIRECTION = "ad";
	public static final String ORIGINAL_CIGAR = SAMTag.OC.name();
	public static final String ORIGINAL_POSITION = SAMTag.OP.name();
	public static final String SPANNING_ASSEMBLY = "sp";
	// Per category aggregations
	public static final String ASSEMBLY_READPAIR_COUNT = "dc";
	public static final String ASSEMBLY_SOFTCLIP_COUNT = "sc";
	public static final String ASSEMBLY_SOFTCLIP_REMOTE_COUNT = "rc";
	public static final String ASSEMBLY_NONSUPPORTING_READPAIR_COUNT = "ec";
	public static final String ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT = "tc";
	public static final String ASSEMBLY_READPAIR_QUAL = "dq";
	public static final String ASSEMBLY_SOFTCLIP_QUAL = "sq";
	public static final String ASSEMBLY_SOFTCLIP_REMOTE_QUAL = "rq";
	public static final String ASSEMBLY_NONSUPPORTING_READPAIR_QUAL = "eq";
	public static final String ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL = "tq";
	// other fields no longer in active use
	public static final String ASSEMBLY_READPAIR_LENGTH_MAX = "dl";
	public static final String ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX = "ds";
	public static final String ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL = "ms";
}
