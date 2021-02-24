package au.edu.wehi.idsv.sam;

public class SamTags {
	/**
	 * Flag indicating whether the "read" is a GRIDSS breakend assembly
	 */
	public static final String IS_ASSEMBLY = "aa";
	/**
	 * Indicates whether this assembly was generated from forward-supporting
	 * or backwards-supporting evidence.
	 *
	 * "b" indicates a backwards assembly with the anchoring bases expected at the end
	 *
	 * "f" indicates a forward assembly with the anchoring bases expected at the start
	 *
	 */
	public static final String ASSEMBLY_DIRECTION = "ad";
	/**
	 * Flag indicating whether the breakend assembly was unanchored.
	 * That is, the assembly did not contain any supporting
	 * soft clipped/split reads and the read pairs were unable to be
	 * assembled all the way back to reference-supporting bases.
	 */
	public static final String UNANCHORED = "ua";
	/**
	 * Strand bias of alignment of supporting soft-clip/split reads.
	 *
	 * Bias is between 0 and 1.
	 * 1 indicates that all assembled SC/SR reads were aligned to the positive strand.
	 * 0 indicates that all assembled SC/SR reads were aligned to the negative strand.
	 * The default value is 0.5
	 */
	public static final String ASSEMBLY_STRAND_BIAS = "sb";
	/**
	 * 0-based start offset of supporting evidence.
	 *
	 * One value per assembled supporting evidence (SC/SR/OEA/DP)
	 * These fields are used to pro-rata assembly support.
	 * This prevents a somatic SV adjacent to a germline indel adjacent being flagged as germline.
	 */
	public static final String ASSEMBLY_EVIDENCE_OFFSET_START = "os";
	/**
	 * 0-based assembly contig start offset of supporting evidence.
	 *
	 * One value per assembled supporting evidence (SC/SR/OEA/DP)
	 * These fields are used to pro-rata assembly support.
	 * This prevents a somatic SV adjacent to a germline indel adjacent being flagged as germline.
	 */
	public static final String ASSEMBLY_EVIDENCE_OFFSET_END = "oe";
	/**
	 * Type of assembly evidence support. See {@link au.edu.wehi.idsv.AssemblyEvidenceSupport.SupportType}
	 *
	 * One value per assembled supporting evidence (SC/SR/OEA/DP)
	 * These fields are used to pro-rata assembly support.
	 * This prevents a somatic SV adjacent to a germline indel adjacent being flagged as germline.
	 */
	public static final String ASSEMBLY_EVIDENCE_TYPE = "et";
	/**
	 * Category of assembly evidence support.
	 * Categories correspond to the 0-based ordinal of the VCF FORMAT support breakdown.
	 * By default there is a 1-1 correspondence between input bams and categories.
	 *
	 * One value per assembled supporting evidence (SC/SR/OEA/DP)
	 * These fields are used to pro-rata assembly support.
	 * This prevents a somatic SV adjacent to a germline indel adjacent being flagged as germline.
	 */
	public static final String ASSEMBLY_EVIDENCE_CATEGORY = "ec";
	/**
	 * QUAL score of assembly evidence support.
	 *
	 * One value per assembled supporting evidence (SC/SR/OEA/DP)
	 * These fields are used to pro-rata assembly support.
	 * This prevents a somatic SV adjacent to a germline indel adjacent being flagged as germline.
	 */
	public static final String ASSEMBLY_EVIDENCE_QUAL = "eq";
	/**
	 * EvidenceID of assembly evidence support. See {@link au.edu.wehi.idsv.EvidenceIdentifierGenerator}
	 *
	 * One value per assembled supporting evidence (SC/SR/OEA/DP)
	 * These fields are used to pro-rata assembly support.
	 * This prevents a somatic SV adjacent to a germline indel adjacent being flagged as germline.
	 */
	public static final String ASSEMBLY_EVIDENCE_EVIDENCEID = "ez";
	/**
	 * Fragment identifier of supporting evidence.
	 * By default, the fragment identifier is just the read name.
	 *
	 * This field is used by {@link au.edu.wehi.idsv.vcf.VcfInfoAttributes#BREAKPOINT_VARIANT_FRAGMENTSVARIANT_FRAGMENTS }
	 * and {@link au.edu.wehi.idsv.vcf.VcfInfoAttributes#BREAKEND_VARIANT_FRAGMENTS } to correctly count the number
	 * of distinct fragments supporting a variant.
	 *
	 * One value per assembled supporting evidence (SC/SR/OEA/DP)
	 * These fields are used to pro-rata assembly support.
	 * This prevents a somatic SV adjacent to a germline indel adjacent being flagged as germline.
	 */
	public static final String ASSEMBLY_EVIDENCE_FRAGMENTID = "ef";
	/**
	 * Maximum length of reads contributing to assembly.
	 */
	public static final String ASSEMBLY_MAX_READ_LENGTH = "dl";
	/***
	 * Pair of integers indicating how many bases at the start/end of the contig were truncated
	 * to ensure the contig anchor did not overrun the contig bounds
	 */
	public static final String ASSEMBLY_ANCHOR_TRUNCATION = "at";
}
