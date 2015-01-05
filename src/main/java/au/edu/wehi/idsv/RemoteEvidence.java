package au.edu.wehi.idsv;

/**
 * Indicates that the underlying source evidence was located from the remote breakend.
 * 
 * @author Daniel Cameron
 *
 */
public interface RemoteEvidence extends DirectedBreakpoint {
	/**
	 * Evidence according to the source evidence at the remote breakend. 
	 * @return evidence from the source remote breakend perspective.
	 */
	DirectedBreakpoint asLocal();
}
