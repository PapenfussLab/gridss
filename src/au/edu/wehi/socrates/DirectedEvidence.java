package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

public interface DirectedEvidence {
	BreakpointDirection getBreakpointDirection();
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
	 * Index of BAM contig this breakpoint is located on.
	 * @return reference index of evidence
	 */
	int getReferenceIndex();
	/**
	 * Start of position of breakpoint location
	 * @return position of breakpoint
	 */
	long getWindowStart();
	long getWindowEnd();
}
