package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

/**
 * Breakpoint evidence including at least one base of sequence at the breakpoint location 
 * @author Daniel Cameron
 *
 */
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
}
