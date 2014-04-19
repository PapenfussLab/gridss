package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

/**
 * Constructs putative directed breakpoints from read evidence
 * providing support for structural variation
 * @author Daniel Cameron
 *
 */
public interface SequentialDirectedBreakpointEvidenceProcessor {
	/**
	 * Gets the evidence for the given position
	 * <p>The position called is guaranteed to be strictly increasing</p>
	 * @param position position to call
	 * @return Evidence for directed breakpoint at the given position. 
	 */
	Iterable<DirectedBreakpoint> getEvidenceAtPosition(int position);
	void addSoftclipSupportingForwardBreakpoint(SAMRecord record);
	void addSoftclipSupportingBackwardBreakpoint(SAMRecord record);
	void addReadPairSupportingForwardBreakpoint(NonReferenceReadPair pair);
	void addReadPairSupportingBackwardBreakpoint(NonReferenceReadPair pair);
}
