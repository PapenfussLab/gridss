package au.edu.wehi.socrates;

import au.edu.wehi.socrates.graph.TrapezoidGraphNode;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

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
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((direction == null) ? 0 : direction.hashCode());
		result = prime * result + end;
		long temp;
		temp = Double.doubleToLongBits(qual);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + referenceIndex;
		result = prime * result + start;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		BreakpointLocation other = (BreakpointLocation) obj;
		if (direction != other.direction)
			return false;
		if (end != other.end)
			return false;
		if (Double.doubleToLongBits(qual) != Double
				.doubleToLongBits(other.qual))
			return false;
		if (referenceIndex != other.referenceIndex)
			return false;
		if (start != other.start)
			return false;
		return true;
	}
	public static Ordering<BreakpointLocation> ByStartEnd = new Ordering<BreakpointLocation>() {
		public int compare(BreakpointLocation o1, BreakpointLocation o2) {
			  return ComparisonChain.start()
			        .compare(o1.referenceIndex, o2.referenceIndex)
			        .compare(o1.start, o2.start)
			        .compare(o1.end, o2.end)
			        .result();
		  }
	};
}
