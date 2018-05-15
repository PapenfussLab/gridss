package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public enum VcfInfoAttributes {
	REFERENCE_READ_COUNT ("REF", 1, VCFHeaderLineType.Integer, "Count of reads mapping across this breakend"),
	REFERENCE_READPAIR_COUNT ("REFPAIR", 1, VCFHeaderLineType.Integer, "Count of reference read pairs spanning this breakpoint supporting the reference allele"),
	REFERENCE_FRAGMENTS ("RF", 1, VCFHeaderLineType.Integer, "Reference fragments. Count of fragments suspporting the reference allele and not the variant allele."),
	CALLED_QUAL ("CQ", 1, VCFHeaderLineType.Float, "Breakpoint quality score before evidence reallocation"),
	BREAKEND_QUAL ("BQ", 1, VCFHeaderLineType.Float, "Quality score of breakend evidence"),
	
	BREAKPOINT_ASSEMBLY_COUNT("AS", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting breakpoint"),
	BREAKPOINT_READPAIR_COUNT("RP", 1, VCFHeaderLineType.Integer, "Count of read pairs supporting breakpoint"),
	BREAKPOINT_SPLITREAD_COUNT("SR", 1, VCFHeaderLineType.Integer, "Count of split reads supporting breakpoint"),
	BREAKPOINT_INDEL_COUNT("IC", 1, VCFHeaderLineType.Integer, "Count of read indels supporting breakpoint"),
	BREAKPOINT_ASSEMBLY_COUNT_REMOTE("RAS", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting breakpoint from remote breakend"),
	BREAKPOINT_ASSEMBLY_COUNT_COMPOUND("CAS", 1, VCFHeaderLineType.Integer, "Count of complex compound breakpoint assemblies supporting breakpoint from elsewhere"),
	
	BREAKPOINT_ASSEMBLY_READPAIR_COUNT("ASRP", 1, VCFHeaderLineType.Integer, "Count of read pairs incorporated into any breakpoint assembly"),
	BREAKPOINT_ASSEMBLY_READ_COUNT("ASSR", 1, VCFHeaderLineType.Integer, "Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies"),
	//BREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT("ASCRP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of read pairs not directly supporting breakpoint incorporated into any breakpoint assembly"),
	//BREAKPOINT_ASSEMBLY_CONSCRIPTED_READ_COUNT("ASCSR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of split, soft clipped or indel-containing reads not directly supporting breakpoint incorporated into any breakpoint assemblies"),
	
	BREAKPOINT_ASSEMBLY_QUAL("ASQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting breakpoint"),
	BREAKPOINT_ASSEMBLY_QUAL_REMOTE("RASQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting breakpoint from remote breakend"),
	BREAKPOINT_ASSEMBLY_QUAL_COMPOUND("CASQ", 1, VCFHeaderLineType.Float, "Quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere"),
	BREAKPOINT_READPAIR_QUAL("RPQ", 1, VCFHeaderLineType.Float, "Quality score of read pairs supporting breakpoint"),
	BREAKPOINT_SPLITREAD_QUAL("SRQ", 1, VCFHeaderLineType.Float, "Quality score of split reads supporting breakpoint"),
	BREAKPOINT_INDEL_QUAL("IQ", 1, VCFHeaderLineType.Float, "Quality score of read indels supporting breakpoint"),
	BREAKPOINT_VARIANT_FRAGMENTS ("VF", 1, VCFHeaderLineType.Integer, "Count of fragments suspporting the variant breakpoint allele and not the reference allele."),

	BREAKEND_ASSEMBLY_COUNT("BA", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting just local breakend"),
	BREAKEND_UNMAPPEDMATE_COUNT("BUM", 1, VCFHeaderLineType.Integer, "Count of read pairs (with one read unmapped) supporting just local breakend"),
	BREAKEND_SOFTCLIP_COUNT("BSC", 1, VCFHeaderLineType.Integer, "Count of soft clips supporting just local breakend"),
	
	BREAKEND_ASSEMBLY_READPAIR_COUNT("BASRP", 1, VCFHeaderLineType.Integer, "Count of read pairs incorporated into any breakend assembly"),
	BREAKEND_ASSEMBLY_READ_COUNT("BASSR", 1, VCFHeaderLineType.Integer, "Count of split, soft clipped or indel-containing reads incorporated into any breakend assemblies"),
	
	BREAKEND_ASSEMBLY_QUAL("BAQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting just local breakend"),
	BREAKEND_UNMAPPEDMATE_QUAL("BUMQ", 1, VCFHeaderLineType.Float, "Quality score of read pairs (with one read unmapped) supporting just local breakend"),
	BREAKEND_SOFTCLIP_QUAL("BSCQ", 1, VCFHeaderLineType.Float, "Quality score of soft clips supporting just local breakend"),
	BREAKEND_VARIANT_FRAGMENTS ("BVF", 1, VCFHeaderLineType.Integer, "Count of fragments providing breakend for the variant allele."),

	CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY ("CIRPOS", 2, VCFHeaderLineType.Integer, "Confidence interval around remote breakend POS for imprecise variants"),
	STRAND_BIAS ("SB", 1, VCFHeaderLineType.Float, "Strand bias of reads supporting the variant."
			+ "1 indicates that breakend bases would be aligned to the positive strand if the reference was changed to the variant allele. "
			+ "0 indicates that breakend bases would be aligned to the negative strand if the reference was changed to the variant allele. "
			+ "For read alignments, this corresponds the strand of the read alignment. "
			+ "A SR that is aligned to the positive strand will have a strand bias of 1. "
			+ "For read pairs in FR orientation, this is the opposite strand to the local anchoring read alignment. "
			+ "A RP that is aligned to the positive strand will have a strand bias of 0 because the mate is providing the support and is expected to align to the negative strand. "
			+ "Note that reads both directly supporting the variant, and supporting via assembly will be double-counted"),
	SELF_INTERSECTING ("SELF", 1, VCFHeaderLineType.Flag, "Indicates a breakpoint is self-intersecting"),
	SUPPORT_INTERVAL ("SI", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Support interval offsets from breakend position in which at least one supporting read/read pair/assembly is mapped."),
	REMOTE_SUPPORT_INTERVAL ("RSI", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Support interval offsets of partner breakend."),
	INEXACT_HOMPOS ("IHOMPOS", 2, VCFHeaderLineType.Integer, "Position of inexact homology"),
	SUPPORT_CIGAR ("SC", 1, VCFHeaderLineType.String, "CIGAR for displaying anchoring alignment of any contributing evidence and microhomologies."),
	BREAKEND_ASSEMBLY_ID ("BEID", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Identifiers of assemblies supporting the variant."),
	BREAKEND_ALIGNMENTS ("BEALN", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Potential alignment locations of breakend sequence in the format chr:start|strand|cigar|mapq. Depending on the alignment information available, strand and mapq may be empty.");
	private final VCFInfoHeaderLine header;
	private final String tag;
	VcfInfoAttributes(String name, String samTag, int count, VCFHeaderLineType type, String description) {
		this(new VCFInfoHeaderLine(name, count, type, description), samTag);
	}
	VcfInfoAttributes(String name, String samTag, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
		this(new VCFInfoHeaderLine(name, count, type, description), samTag);
	}
	VcfInfoAttributes(String name, int count, VCFHeaderLineType type, String description) {
		this(name, null, count, type, description);
	}
	VcfInfoAttributes(String name, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
		this(name, null, count, type, description);
	}
	VcfInfoAttributes(VCFInfoHeaderLine header, String samTag) {
		this.header = header;
		this.tag = samTag;
	}
	public VCFInfoHeaderLine infoHeader() { return header; }
	public String attribute() { return header != null ? header.getID() : null; }
	/**
	 * Gets the attribute for the given key
	 * @param key VCF info field name
	 * @return corresponding attribute, null if no idsv-specific attribute with the given name is defined
	 */
	public static VcfInfoAttributes getAttributefromKey(String key) {
		for (VcfInfoAttributes a : values()) {
			if (a.attribute().equals(key)) return a;
		}
		return null;
	}
	/**
	 * SAM tag for attribute when persisted to BAM
	 * @return
	 */
	public String samTag() {
		return tag;
	}
}
