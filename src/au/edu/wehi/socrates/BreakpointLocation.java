package au.edu.wehi.socrates;

/**
 * Positional location of a breakpoint that is consistent with the given evidence
 * @author Daniel Cameron
 *
 */
public class BreakpointLocation {
	/**
	 * First possible position of breakpoint is immediately after this 1-based genomic coordinate
	 */
	public final long start;
	/**
	 * Last possible position of breakpoint is immediately after this 1-based genomic coordinate
	 */
	public final long end;
	/**
	 * Breakpoint position is on this contig
	 */
	public final int referenceIndex;
	/**
	 * Breakpoint is in the given direction
	 */
	public final BreakpointDirection direction;
	/**
	 * Phred-like quality score of this evidence
	 */
	public final float qual;
	public BreakpointLocation(int referenceIndex, BreakpointDirection direction, long start, long end, float qual) {
		this.referenceIndex = referenceIndex;
		this.direction = direction;
		this.start = start;
		this.end = end;
		this.qual = qual;
	}
}
