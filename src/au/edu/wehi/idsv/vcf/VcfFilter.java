package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFFilterHeaderLine;

/**
 * VCF filters used by idsv
 * @author cameron.d
 *
 */
public enum VcfFilter {
	ASSEMBLY_TOO_SHORT ("ASSEMBLY_TOO_SHORT", "This assembly is shorter than a read length"),
	ASSEMBLY_SINGLE_READ ("ASSEMBLY_SINGLE_READ", "Only a single read contributes to this assembly."),
	ASSEMBLY_REF ("ASSEMBLY_REF", "Assembly corresponds to reference allele");

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
}
