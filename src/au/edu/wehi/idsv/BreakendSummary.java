package au.edu.wehi.idsv;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

/**
 * Positional location of a breakpoint that is consistent with the given evidence
 * @author Daniel Cameron
 *
 */
public class BreakendSummary {
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
	public final BreakendDirection direction;
	/**
	 * Phred-like quality score of this evidence
	 */
	public final EvidenceMetrics evidence;
	public BreakendSummary(int referenceIndex, BreakendDirection direction, int start, int end, EvidenceMetrics evidence) {
		if (referenceIndex < 0) {
			throw new IllegalArgumentException("Reference index must be valid");
		}
		this.referenceIndex = referenceIndex;
		this.direction = direction;
		this.start = start;
		this.end = end;
		this.evidence = evidence;
	}
	public BreakendSummary(BreakendSummary location, EvidenceMetrics evidence) {
		this.referenceIndex = location.referenceIndex;
		this.direction = location.direction;
		this.start = location.start;
		this.end = location.end;
		this.evidence = evidence;
	}
	protected static String toString(int referenceIndex, int start, int end) {
		return String.format("%d:%d-%d", referenceIndex, start, end);
	}
	protected static String toString(BreakendDirection direction, int referenceIndex, int start, int end) {
		if (direction == BreakendDirection.Forward) return toString(referenceIndex, start, end) + ">";
		return "<" + toString(referenceIndex, start, end);
	}
	@Override
	public String toString() {
		return String.format("%s %s", toString(direction, referenceIndex, start, end), evidence);
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((direction == null) ? 0 : direction.hashCode());
		result = prime * result + end;
		long temp;
		temp = evidence.hashCode();
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
		BreakendSummary other = (BreakendSummary) obj;
		if (direction != other.direction)
			return false;
		if (end != other.end)
			return false;
		if (!evidence.equals(other.evidence))
			return false;
		if (referenceIndex != other.referenceIndex)
			return false;
		if (start != other.start)
			return false;
		return true;
	}
	public static Ordering<BreakendSummary> ByStartEnd = new Ordering<BreakendSummary>() {
		public int compare(BreakendSummary o1, BreakendSummary o2) {
			  return ComparisonChain.start()
			        .compare(o1.referenceIndex, o2.referenceIndex)
			        .compare(o1.start, o2.start)
			        .compare(o1.end, o2.end)
			        .result();
		  }
	};
	public static Ordering<BreakendSummary> ByEndStart = new Ordering<BreakendSummary>() {
		public int compare(BreakendSummary o1, BreakendSummary o2) {
			  return ComparisonChain.start()
			        .compare(o1.referenceIndex, o2.referenceIndex)
			        .compare(o1.end, o2.end)
			        .compare(o1.start, o2.start)
			        .result();
		  }
	};
	/**
	 * This breakend is fully contained by the given breakend
	 * @param other
	 * @return
	 */
	public boolean containedBy(BreakendSummary other) {
		return this.referenceIndex == other.referenceIndex &&
				this.direction == other.direction &&
				this.start >= other.start &&
				this.end <= other.end;
	}
	/**
	 * This breakend shares at least one position with the given breakend
	 * @param loc
	 * @return
	 */
	public boolean overlaps(BreakendSummary loc) {
		return breakendOverlaps(loc);
	}
	protected boolean breakendOverlaps(BreakendSummary loc) {
		return this.referenceIndex == loc.referenceIndex &&
				this.direction == loc.direction &&
				((this.start <= loc.start && this.end >= loc.start) ||
				 (this.start >= loc.start && this.start <= loc.end));
	}
}
