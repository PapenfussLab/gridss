package au.edu.wehi.socrates;

import org.broadinstitute.variant.variantcontext.VariantContext;

public interface DirectedBreakpoint {
	/**
	 * Unique breakpoint identifier.
	 * This identifier is used as the FASTQ sequence identifier
	 * when performing realignment and must be a valid FASTQ
	 * sequence identifier as well as being unique across all
	 * evidence processors.</p>
	 * @return Unique breakpoint identifier string
	 */
	String getBreakpointID();
	/**
	 * Index of BAM contig this breakpoint is located on.
	 * @return reference index of breakpoint
	 */
	String getChr();
	/**
	 * Start of position of breakpoint location
	 * @return position of breakpoint
	 */
	int getStart();
	int getEnd();
	BreakpointDirection getBreakpointDirection();
	byte[] getBreakpointSequence();
	byte[] getBreakpointQuality();
	VariantContext getVariantContext();
}
