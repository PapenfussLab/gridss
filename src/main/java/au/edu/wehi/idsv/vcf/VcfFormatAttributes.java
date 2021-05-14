package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

public enum VcfFormatAttributes {
	REFERENCE_READ_COUNT ("REF", 1, VCFHeaderLineType.Integer, "Count of reads mapping across this breakend"),
	REFERENCE_READPAIR_COUNT ("REFPAIR", 1, VCFHeaderLineType.Integer, "Count of reference read pairs spanning this breakend supporting the reference allele"),
	//REFERENCE_FRAGMENTS ("RF", 1, VCFHeaderLineType.Integer, "Reference fragments. Count of fragments supporting the reference allele and not the variant allele."),
	BREAKEND_QUAL ("BQ", 1, VCFHeaderLineType.Float, "Quality score of breakend evidence after evidence reallocation"),
	BREAKPOINT_QUAL ("QUAL", 1, VCFHeaderLineType.Float, "Quality score of breakend evidence after evidence reallocation"),
	
	BREAKPOINT_READPAIR_COUNT("RP", 1, VCFHeaderLineType.Integer, "Count of read pairs supporting breakpoint per category"),
	BREAKPOINT_SPLITREAD_COUNT("SR", 1, VCFHeaderLineType.Integer, "Count of split reads supporting breakpoint per category"),
	BREAKPOINT_INDEL_COUNT("IC", 1, VCFHeaderLineType.Integer, "Count of read indels supporting breakpoint per category"),
	BREAKPOINT_ASSEMBLY_READPAIR_COUNT("ASRP", 1, VCFHeaderLineType.Integer, "Count of read pairs incorporated into any breakpoint assembly"),
	BREAKPOINT_ASSEMBLY_READ_COUNT("ASSR", 1, VCFHeaderLineType.Integer, "Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies"),
	BREAKPOINT_ASSEMBLY_QUAL("ASQ", 1, VCFHeaderLineType.Float, "Pro-rata quality score contribution of assemblies supporting breakpoint"),
	BREAKPOINT_ASSEMBLY_QUAL_REMOTE("RASQ", 1, VCFHeaderLineType.Float, "Pro-rata quality score contribution of assemblies supporting breakpoint from remote breakend"),
	BREAKPOINT_ASSEMBLY_QUAL_COMPOUND("CASQ", 1, VCFHeaderLineType.Float, "Pro-rata quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere"),
	BREAKPOINT_READPAIR_QUAL("RPQ", 1, VCFHeaderLineType.Float, "Quality score of read pairs supporting breakpoint per category"),
	BREAKPOINT_SPLITREAD_QUAL("SRQ", 1, VCFHeaderLineType.Float, "Quality score of split reads supporting breakpoint per category"),
	BREAKPOINT_INDEL_QUAL("IQ", 1, VCFHeaderLineType.Float, "Quality score of read indels supporting breakpoint per category"),
	BREAKPOINT_VARIANT_FRAGMENTS ("VF", 1, VCFHeaderLineType.Integer, "Count of fragments supporting the variant breakpoint allele and not the reference allele."),
	BREAKPOINT_ASSEMBLED_NONSUPPORTING_READPAIR_COUNT("BANRP", 1, VCFHeaderLineType.Integer, "Count of read pairs at this breakend assembled into a contig that does not support the breakpoint."),
	BREAKPOINT_ASSEMBLED_NONSUPPORTING_SPLITREAD_COUNT("BANSR", 1, VCFHeaderLineType.Integer, "Count of split reads at this breakend assembled into a contig that does not support the breakpoint."),
	BREAKPOINT_ASSEMBLED_NONSUPPORTING_READPAIR_QUAL("BANRPQ", 1, VCFHeaderLineType.Float, "Quality score of read pairs at this breakend assembled into a contig that does not support the breakpoint."),
	BREAKPOINT_ASSEMBLED_NONSUPPORTING_SPLITREAD_QUAL("BANSRQ", 1, VCFHeaderLineType.Float, "Quality score of split reads at this breakend assembled into a contig that does not support the breakpoint."),

	BREAKEND_UNMAPPEDMATE_COUNT("BUM", 1, VCFHeaderLineType.Integer, "Count of read pairs (with one read unmapped) supporting just local breakend per category"),
	BREAKEND_SOFTCLIP_COUNT("BSC", 1, VCFHeaderLineType.Integer, "Count of soft clips supporting just local breakend per category"),
	BREAKEND_ASSEMBLY_READPAIR_COUNT("BASRP", 1, VCFHeaderLineType.Integer, "Count of read pairs incorporated into any breakend assembly"),
	BREAKEND_ASSEMBLY_READ_COUNT("BASSR", 1, VCFHeaderLineType.Integer, "Count of split, soft clipped or indel-containing reads incorporated into any breakend assemblies"),
	BREAKEND_ASSEMBLY_QUAL("BAQ", 1, VCFHeaderLineType.Float, "Pro-rata quality score contribution of assemblies supporting just local breakend"),
	BREAKEND_UNMAPPEDMATE_QUAL("BUMQ", 1, VCFHeaderLineType.Float, "Quality score of read pairs (with one read unmapped) supporting just local breakend per category"),
	BREAKEND_SOFTCLIP_QUAL("BSCQ", 1, VCFHeaderLineType.Float, "Quality score of soft clips supporting just local breakend per category"),
	BREAKEND_VARIANT_FRAGMENTS ("BVF", 1, VCFHeaderLineType.Integer, "Count of fragments providing breakend for the variant allele.");

	private final VCFFormatHeaderLine header;
	private final String tag;
	VcfFormatAttributes(String name, String samTag, int count, VCFHeaderLineType type, String description) {
		this(new VCFFormatHeaderLine(name, count, type, description), samTag);
	}
	VcfFormatAttributes(String name, String samTag, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
		this(new VCFFormatHeaderLine(name, count, type, description), samTag);
	}
	VcfFormatAttributes(String name, int count, VCFHeaderLineType type, String description) {
		this(name, null, count, type, description);
	}
	VcfFormatAttributes(String name, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
		this(name, null, count, type, description);
	}
	VcfFormatAttributes(VCFFormatHeaderLine header, String samTag) {
		this.header = header;
		this.tag = samTag;
	}
	public VCFFormatHeaderLine formatHeader() { return header; }
	public String attribute() { return header != null ? header.getID() : null; }
	/**
	 * Gets the attribute for the given key
	 * @param key VCF info field name
	 * @return corresponding attribute, null if no idsv-specific attribute with the given name is defined
	 */
	public static VcfFormatAttributes getAttributefromKey(String key) {
		for (VcfFormatAttributes a : values()) {
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
