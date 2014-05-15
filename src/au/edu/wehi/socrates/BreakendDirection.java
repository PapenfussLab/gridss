package au.edu.wehi.socrates;

/**
 * Direction of breakpoint on the given chromosome.
 * <p>The positional relationship between the direct and breakpoint position
 * matches that used by VCF.</p>
 * @author Daniel Cameron
 *
 */
public enum BreakendDirection {
	/**
	 * The breakpoint includes bases at and before the breakpoint position
	 * but not after.
	 * 
	 * AAAA]?]
	 */
	Forward,
	/**
	 * The breakpoint iincludes bases after but not at the breakpoint position.
	 * 
	 * N[?[AAAA
	 */
	Backward,
}
