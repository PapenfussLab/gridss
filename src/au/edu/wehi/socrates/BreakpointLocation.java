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
	protected static String toString(int referenceIndex, int start, int end) {
		return String.format("%d:%d-%d");
	}
	protected static String toString(BreakpointDirection direction, int referenceIndex, int start, int end) {
		if (direction == BreakpointDirection.Forward) return toString(referenceIndex, start, end) + ">";
		return "<" + toString(referenceIndex, start, end);
	}
	@Override
	public String toString() {
		return String.format("%s %.1f", toString(direction, referenceIndex, start, end), qual);
	}
}
