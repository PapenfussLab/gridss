package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

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
