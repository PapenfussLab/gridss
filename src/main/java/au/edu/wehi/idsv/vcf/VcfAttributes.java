package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public enum VcfAttributes {
	REFERENCE_READ_COUNT ("REF", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of reads mapping across this breakend per category"),
	REFERENCE_READPAIR_COUNT ("REFPAIR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of reference read pairs spanning this breakpoint supporting the reference allele per category"),
	//SOMATIC_P_VALUE ("SPV", 1, VCFHeaderLineType.Float, "Somatic p-value. Assume first category is normal, second is tumour"),
	CALLED_QUAL ("CQ", 1, VCFHeaderLineType.Float, "Breakpoint quality score before evidence reallocation"),
	//BREAKPOINT_QUAL ("BPQUAL", 1, VCFHeaderLineType.Float, "Quality score of breakpoint evidence after evidence reallocation"), // QUAL field
	BREAKEND_QUAL ("BQ", 1, VCFHeaderLineType.Float, "Quality score of breakend evidence after evidence reallocation"),
	
	BREAKPOINT_ASSEMBLY_COUNT("AS", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting breakpoint"),
	BREAKPOINT_READPAIR_COUNT("RP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of read pairs supporting breakpoint per category"),
	BREAKPOINT_SPLITREAD_COUNT("SR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of split reads supporting breakpoint per category"),
	BREAKPOINT_ASSEMBLY_COUNT_REMOTE("RAS", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting breakpoint from remote breakend"),
	BREAKPOINT_SPLITREAD_COUNT_REMOTE("RSR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of split reads supporting breakpoint from remote breakend per category"),
	
	BREAKPOINT_ASSEMBLY_READPAIR_COUNT("ASRP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of read pairs incorporated into any breakpoint assembly"),
	BREAKPOINT_ASSEMBLY_SPLITREAD_COUNT("ASSR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of split or soft clipped reads incorporated into any breakpoint assemblies"),
	BREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT("ASCRP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of read pairs not directly supporting breakpoint incorporated into any breakpoint assembly"),
	BREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT("ASCSR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of split or soft clipped reads not directly supporting breakpoint incorporated into any breakpoint assemblies"),
	
	BREAKPOINT_ASSEMBLY_QUAL("ASQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting breakpoint"),
	BREAKPOINT_READPAIR_QUAL("RPQ", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Quality score of read pairs supporting breakpoint per category"),
	BREAKPOINT_SPLITREAD_QUAL("SRQ", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Quality score of split reads supporting breakpoint per category"),
	BREAKPOINT_ASSEMBLY_QUAL_REMOTE("RASQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting breakpoint from remote breakend"),
	BREAKPOINT_SPLITREAD_QUAL_REMOTE("RSRQ", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Quality score of split reads supporting breakpoint from remote breakend per category"),

	BREAKEND_ASSEMBLY_COUNT("BA", 1, VCFHeaderLineType.Integer, "Count of assemblies supporting just local breakend"),
	BREAKEND_UNMAPPEDMATE_COUNT("BUM", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of read pairs supporting just local breakend per category"),
	BREAKEND_SOFTCLIP_COUNT("BSC", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Count of soft clips supporting just local breakend per category"),
	BREAKEND_ASSEMBLY_QUAL("BAQ", 1, VCFHeaderLineType.Float, "Quality score of assemblies supporting just local breakend"),
	BREAKEND_UNMAPPEDMATE_QUAL("BUMQ", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Quality score of read pairs supporting just local breakend per category"),
	BREAKEND_SOFTCLIP_QUAL("BSCQ", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Quality score of soft clips supporting just local breakend per category"),

	CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY ("CIRPOS", 2, VCFHeaderLineType.Integer, "Confidence interval around remote breakend POS for imprecise variants"),
	
	SELF_INTERSECTING ("SELF", 1, VCFHeaderLineType.Flag, "Indicates a breakpoint is self-intersecting"),
	ANCHOR_CIGAR ("CIGAR", 1, VCFHeaderLineType.String, "CIGAR for displaying anchoring alignment of any contributing evidence and microhomologies."),
	TEST ("TEST", 1, VCFHeaderLineType.String, "Placeholder field used for regression testing.");
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
