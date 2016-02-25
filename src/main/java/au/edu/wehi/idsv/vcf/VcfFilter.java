package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFFilterHeaderLine;

/**
 * VCF filters used by idsv
 * @author Daniel Cameron
 *
 */
public enum VcfFilter {
	ASSEMBLY_TOO_SHORT ("ASSEMBLY_TOO_SHORT", "This assembly is shorter than a read length"),
	ASSEMBLY_TOO_FEW_READ ("ASSEMBLY_TOO_FEW_READ", "Not enough reads contribute to this assembly as specified by 'assembly.minReads'"),
	REFERENCE_ALLELE ("REF", "Breakpoint corresponds to reference allele"),
	SMALL_EVENT ("SMALL_EVENT", "Event size is smaller than the minimum reportable size specified by 'variantcalling.minSize'"),
	LOW_BREAKPOINT_SUPPORT ("LOW_BREAKPOINT_SUPPORT", "Does not reach the required threshold quality for calling as specified by 'variantcalling.minScore'"),
	SINGLE_SUPPORT ("SINGLE_SUPPORT", "Supported by a single read or read pair only."),
	//ASSEMBLY_REMOTE ("ASSEMBLY_REMOTE", "All support for the given breakpoint comes from elsewhere."),
	NO_ASSEMBLY ("NO_ASSEMBLY", "No assembly supporting this variant could be found."),
	SINGLE_ASSEMBLY ("SINGLE_ASSEMBLY", "Only one side of the breakpoint could be assembled."),
	ASSEMBLY_ONLY ("ASSEMBLY_ONLY", "Variant is supported only by assembly evidence."),
	LOW_QUAL ("LOW_QUAL", "Low quality call as specified by 'variantcalling.lowQuality'");

    private final VCFFilterHeaderLine filterheader;
	VcfFilter(String name, String description) {
		this.filterheader = new VCFFilterHeaderLine(name, description);
	}
	/**
	 * Name of the VCF filter
	 * @return
	 */
	public String filter() { return filterheader.getID(); }
	/**
	 * VCF metadata header for this filter 
	 * @return
	 */
	public VCFFilterHeaderLine header() { return filterheader; }	
	public static VcfFilter get(String filter) {
		for (VcfFilter f : values()) {
			if (f.filter().equals(filter)) {
				return f;
			}
		}
		throw new IllegalArgumentException(String.format("Filter %s not found.", filter));
	}
}
