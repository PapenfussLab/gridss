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
	 * EvidenceID of uncategorised assembly components
	 */
	public static final String EVIDENCEID_UNCATEGORISED = "ez";
	/**
	 * EvidenceID of constitute assembly evidence
	 */
	public static final String EVIDENCEID_ASSEMBLY = "ea";
	/**
	 * EvidenceID of constitute remote assembly evidence
	 */
	public static final String EVIDENCEID_ASSEMBLY_REMOTE = "er";
	/**
	 * EvidenceID of constitute breakend assembly evidence
	 */
	public static final String EVIDENCEID_BREAKEND_ASSEMBLY = "eb";
	/**
	 * EvidenceID of constitute discordant pair evidence
	 */
	public static final String EVIDENCEID_PREFIX_DISCORDANT_PAIR = "d";
	/**
	 * EvidenceID of constitute unmapped mate evidence
	 */
	public static final String EVIDENCEID_PREFIX_UNMAPPED_MATE = "u";
	/**
	 * EvidenceID of constitute split read evidence
	 */
	public static final String EVIDENCEID_PREFIX_SPLIT_READ = "s";
	/**
	 * EvidenceID of constitute split read remote evidence
	 */
	public static final String EVIDENCEID_PREFIX_SPLIT_READ_REMOTE = "t";
	/**
	 * EvidenceID of constitute soft clip evidence
	 */
	public static final String EVIDENCEID_PREFIX_SOFT_CLIP = "c";
	
	public static final String ASSEMBLY_READPAIR_COUNT = "dp";
	public static final String ASSEMBLY_SOFTCLIP_COUNT = "sc";
	public static final String ASSEMBLY_REMOTE_COUNT = "rc";
	public static final String ASSEMBLY_NONSUPPORTING_COUNT = "ns";
	public static final String ASSEMBLY_READPAIR_LENGTH_MAX = "rl";
	public static final String ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX = "cl";
	public static final String ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL = "ct";
	public static final String ASSEMBLY_READPAIR_QUAL = "qp";
	public static final String ASSEMBLY_SOFTCLIP_QUAL = "qs";
	public static final String ASSEMBLY_REMOTE_QUAL = "qr";
	public static final String ASSEMBLY_NONSUPPORTING_QUAL = "qn";
	public static final String ASSEMBLY_DIRECTION = "ad";
	public static final String ORIGINAL_CIGAR = "oc";
	public static final String ORIGINAL_POSITION = "op";
	public static final String SPANNING_ASSEMBLY = "sp";
	
}
