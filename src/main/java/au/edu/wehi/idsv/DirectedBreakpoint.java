package au.edu.wehi.idsv;


/**
 * Breakpoint evidence including at least one base of sequence at the breakpoint location 
 * @author Daniel Cameron
 *
 */
public interface DirectedBreakpoint extends DirectedEvidence {
	public BreakpointSummary getBreakendSummary();
	int getRemoteMapq();
	int getRemoteBaseLength();
	int getRemoteBaseCount();
	int getRemoteMaxBaseQual();
	int getRemoteTotalBaseQual();
}
