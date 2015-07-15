package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFFilterHeaderLine;

/**
 * VCF filters used by idsv
 * @author cameron.d
 *
 */
public enum VcfFilter {
	ASSEMBLY_TOO_SHORT ("ASSEMBLY_TOO_SHORT", "This assembly is shorter than a read length"),
	ASSEMBLY_TOO_FEW_READ ("ASSEMBLY_TOO_FEW_READ", "Not enough reads contribute to this assembly."),
	ASSEMBLY_INCORRECT_DIRECTION ("ASSEMBLY_INCORRECT_DIRECTION", "The breakend direction does not match the breakend direction of the assembled reads."),
	REFERENCE_ALLELE ("REF", "Breakpoint corresponds to reference allele"),
	SMALL_INDEL ("SMALL_INDEL", "Variant is a small indel likely to be a false positive."),
	LOW_BREAKPOINT_SUPPORT ("LOW_BREAKPOINT_SUPPORT", "Insufficent support for the given breakpoint"),
	ASSEMBLY_REMOTE ("ASSEMBLY_REMOTE", "All support for the given breakpoint comes from elsewhere."),
	NO_ASSEMBLY ("NO_ASSEMBLY", "No assembly supporting this variant could be found. This variant is of low confidence.");
	//SINGLE_ASSEMBLY ("SINGLE_ASSEMBLY", "Only one side of the breakpoint could be assembled. This variant is of medium confidence.");

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
