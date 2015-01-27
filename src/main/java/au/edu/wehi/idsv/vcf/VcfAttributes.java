package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public enum VcfAttributes {
	REFERENCE_READ_COUNT ("REF", 2, VCFHeaderLineType.Integer, "Count of reads mapping across this breakend (Normal,Tumour)"),
	REFERENCE_READPAIR_COUNT ("REFPAIR", 2, VCFHeaderLineType.Integer, "Count of reference read pairs spanning this breakpoint supporting the reference allele (Normal,Tumour)"),
	SOMATIC_P_VALUE ("SPV", 1, VCFHeaderLineType.Float, "Somatic p-value"),
	CALLED_QUAL ("CQ", 1, VCFHeaderLineType.Float, "Breakpoint quality score before evidence reallocation"),
	//BREAKPOINT_QUAL ("BPQUAL", 1, VCFHeaderLineType.Float, "Quality score of breakpoint evidence after evidence reallocation"), // QUAL field
	BREAKEND_QUAL ("BQ", 1, VCFHeaderLineType.Float, "Quality score of breakend evidence after evidence reallocation"),
	
	BREAKPOINT_ASSEMBLY_COUNT("AS", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting breakpoint"),
	BREAKPOINT_READPAIR_COUNT("RP", 2, VCFHeaderLineType.Integer, "Count of read pairs supporting breakpoint (Normal,Tumour)"),
	BREAKPOINT_SOFTCLIP_COUNT("SC", 2, VCFHeaderLineType.Integer, "Count of soft clips supporting breakpoint (Normal,Tumour)"),
	BREAKPOINT_ASSEMBLY_COUNT_REMOTE("RAS", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting breakpoint from remote breakend"),
	BREAKPOINT_SOFTCLIP_COUNT_REMOTE("RSC", 2, VCFHeaderLineType.Integer, "Count of soft clips supporting breakpoint from remote breakend (Normal,Tumour)"),
	
	BREAKPOINT_ASSEMBLY_QUAL("ASQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting breakpoint"),
	BREAKPOINT_READPAIR_QUAL("RPQ", 2, VCFHeaderLineType.Float, "Quality score of read pairs supporting breakpoint (Normal,Tumour)"),
	BREAKPOINT_SOFTCLIP_QUAL("SCQ", 2, VCFHeaderLineType.Float, "Quality score of soft clips supporting breakpoint (Normal,Tumour)"),
	BREAKPOINT_ASSEMBLY_QUAL_REMOTE("RASQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting breakpoint from remote breakend"),
	BREAKPOINT_SOFTCLIP_QUAL_REMOTE("RSCQ", 2, VCFHeaderLineType.Float, "Quality score of soft clips supporting breakpoint from remote breakend (Normal,Tumour)"),

	BREAKEND_ASSEMBLY_COUNT("BAS", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting just local breakend"),
	BREAKEND_READPAIR_COUNT("BRP", 2, VCFHeaderLineType.Integer, "Count of read pairs supporting just local breakend (Normal,Tumour)"),
	BREAKEND_SOFTCLIP_COUNT("BSC", 2, VCFHeaderLineType.Integer, "Count of soft clips supporting just local breakend (Normal,Tumour)"),

	BREAKEND_ASSEMBLY_QUAL("BASQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting just local breakend"),
	BREAKEND_READPAIR_QUAL("BRPQ", 2, VCFHeaderLineType.Float, "Quality score of read pairs supporting just local breakend (Normal,Tumour)"),
	BREAKEND_SOFTCLIP_QUAL("BSCQ", 2, VCFHeaderLineType.Float, "Quality score of soft clips supporting just local breakend (Normal,Tumour)"),

	CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY ("CIRPOS", 2, VCFHeaderLineType.Integer, "Confidence interval around remote breakend POS for imprecise variants"),
	
	SELF_INTERSECTING ("SELF", 1, VCFHeaderLineType.Flag, "Indicates a breakpoint is self-intersecting"),
	
	TRUTH_MATCHES ("TRUTH_MATCHES", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Variant ID of all true variants as defined in the supplied truth set"),
	TRUTH_MISREALIGN ("TRUTH_MISREALIGN", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Variant ID of the true variant for all variants where the remote breakend location does not match the true variant as defined in the supplied truth set");
	private final VCFInfoHeaderLine header;
	private final String tag;
	VcfAttributes(String name, String samTag, int count, VCFHeaderLineType type, String description) {
		this(new VCFInfoHeaderLine(name, count, type, description), samTag);
	}
	VcfAttributes(String name, String samTag, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
		this(new VCFInfoHeaderLine(name, count, type, description), samTag);
	}
	VcfAttributes(String name, int count, VCFHeaderLineType type, String description) {
		this(name, null, count, type, description);
	}
	VcfAttributes(String name, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
		this(name, null, count, type, description);
	}
	VcfAttributes(VCFInfoHeaderLine header, String samTag) {
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
	public static VcfAttributes getAttributefromKey(String key) {
		for (VcfAttributes a : values()) {
			if (a.attribute() == key) return a;
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
