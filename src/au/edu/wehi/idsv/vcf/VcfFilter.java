package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFFilterHeaderLine;

public enum VcfFilter {
	ASSEMBLY_TOO_SHORT ("ASSEMBLED_SINGLE_READ", "This assembly is shorter than a read length"),
	ASSEMBLY_SINGLE_READ ("ASSEMBLED_SINGLE_READ", "Only a single read contributes to this assembly."),
	ASSEMBLY_REF ("ASSEMBLED_REF", "Assembly corresponds to reference allele");

    private final VCFFilterHeaderLine filterheader;
	VcfFilter(String name, String description) {
		this.filterheader = new VCFFilterHeaderLine(name, description);
	}
	public VCFFilterHeaderLine header() { return filterheader; }	
}
