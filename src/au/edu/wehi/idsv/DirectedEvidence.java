package au.edu.wehi.idsv;

public interface DirectedEvidence {
	/**
	 * Location of breakpoints consistent with the given evidence.
	 * If the destination of the breakpoint is known, a @see BreakpointSummary
	 * should be returned. 
	 * @return breakpoint locations implied by this evidence
	 */
	BreakendSummary getBreakendSummary();
	/**
	 * Unique breakpoint identifier.
	 * This identifier is used as the FASTQ sequence identifier
	 * when performing realignment and must be a valid FASTQ
	 * sequence identifier as well as being unique across all
	 * evidence processors.</p>
	 * @return Unique breakpoint identifier string
	 */
	String getEvidenceID();
}
