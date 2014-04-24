package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

public interface RealignedDirectedBreakpoint extends DirectedBreakpoint {
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
