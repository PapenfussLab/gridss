package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

public interface DirectedBreakpoint extends DirectedEvidence {
	/**
	 * Gets the breakpoint sequence excluding anchor.
	 * @return breakpoint sequence
	 */
	byte[] getBreakpointSequence();
	/**
	 * Gets the breakpoint sequence quality
	 * @return 0-based phred-like quality scores
	 */
	byte[] getBreakpointQuality();
	/**
	 * Sets the realigned breakpoint consensus sequence
	 * @param realigned
	 */
	void setRealigned(SAMRecord realigned);
	/**
	 * Gets the breakpoint sequence realignment location
	 * @return 
	 */
	SAMRecord getRealigned();
}
