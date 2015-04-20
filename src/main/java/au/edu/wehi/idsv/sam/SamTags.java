package au.edu.wehi.idsv.sam;

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
	public static final String ASSEMBLY_COMPONENT_EVIDENCEID = "es";
	public static final String ASSEMBLY_BASE_COUNT = "bc";
	public static final String ASSEMBLY_READPAIR_COUNT = "dp";
	public static final String ASSEMBLY_SOFTCLIP_COUNT = "sc";
	public static final String ASSEMBLY_REMOTE_COUNT = "rc";
	public static final String ASSEMBLY_NONSUPPORTING_COUNT = "ns";
	public static final String ASSEMBLY_READPAIR_LENGTH_MAX = "rl";
	public static final String ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX = "sl";
	public static final String ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL = "st";
	public static final String ASSEMBLY_READPAIR_QUAL = "qp";
	public static final String ASSEMBLY_SOFTCLIP_QUAL = "qs";
	public static final String ASSEMBLY_REMOTE_QUAL = "qr";
	public static final String ASSEMBLY_NONSUPPORTING_QUAL = "qn";
	public static final String ASSEMBLY_DIRECTION = "ad";
	public static final String ORIGINAL_CIGAR = "oc";
	public static final String ORIGINAL_POSITION = "op";
	public static final String SPANNING_ASSEMBLY = "sp";
	
}
