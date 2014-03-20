package au.edu.wehi.socrates;

public enum BreakpointDirection {
	/**
	 * The breakpoint is an adjacency of reference bases immediately before the breakpoint position
	 * and non-reference sequence after the breakpoint position
	 * 
	 * AAAA]?]
	 */
	Forward,
	/**
	 * The breakpoint is an adjacency of reference bases immediately after the breakpoint position
	 * and non-reference sequence before the breakpoint position
	 * 
	 * [?[AAAA
	 */
	Backward,
}
