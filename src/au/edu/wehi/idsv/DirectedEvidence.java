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
	/**
	 * Source of this evidence
	 * @return Source providing this evidence
	 */
	EvidenceSource getEvidenceSource();
	
	/**
	 * Log-likelihood ratio of existence of a structural variation supporting allele vs all reference alleles
	 * @return Log-likelihood ratio
	 */
	float getPhredLogLikelihoodRatio();
	/**
	 * MAPQ of SV-supporting evidence mapped to the reference 
	 * @return
	 */
	int getLocalMapq();
	/**
	 * Length of reference mapped sequence
	 * @return
	 */
	int getLocalBaseLength();
	/**
	 * Total number of reference mapped bases
	 * Note: this will match {@link getLocalBaseLength()} for raw reads 
	 * @return
	 */
	int getLocalBaseCount();
	/**
	 * Maximum base quality of reference mapped bases
	 * @return
	 */
	int getLocalMaxBaseQual();
	/**
	 * Total base quality of reference mapped bases
	 * @return
	 */
	int getLocalTotalBaseQual();
}
