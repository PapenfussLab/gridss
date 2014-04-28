package au.edu.wehi.socrates;

/**
 * Positional location of a breakpoint that is consistent with the given evidence
 * @author Daniel Cameron
 *
 */
public class BreakpointLocation {
	/**
	 * First possible mapped position adjacent to breakpoint in 1-based genomic coordinate
	 */
	public final int start;
	/**
	 * Last possible mapped position adjacent to breakpoint in 1-based genomic coordinate
	 */
	public final int end;
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
	public final double qual;
	public BreakpointLocation(int referenceIndex, BreakpointDirection direction, int start, int end, double qual) {
		this.referenceIndex = referenceIndex;
		this.direction = direction;
		this.start = start;
		this.end = end;
		this.qual = qual;
	}
	public BreakpointLocation(BreakpointLocation location, double qual) {
		this.referenceIndex = location.referenceIndex;
		this.direction = location.direction;
		this.start = location.start;
		this.end = location.end;
		this.qual = qual;
	}
}
