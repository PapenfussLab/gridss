package au.edu.wehi.idsv;

import java.math.RoundingMode;

import htsjdk.samtools.SAMSequenceDictionary;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.math.IntMath;

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
	public BreakendSummary(int referenceIndex, BreakendDirection direction, int start, int end) {
		if (referenceIndex < 0) {
			throw new IllegalArgumentException("Reference index must be valid");
		}
		if (start <= 0) {
			throw new IllegalArgumentException("start must be valid");
		}
		if (end < start) {
			throw new IllegalArgumentException("end must be at or after start");
		}
		this.referenceIndex = referenceIndex;
		this.direction = direction;
		this.start = start;
		this.end = end;
	}
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
	/**
	 * Returns the overlap of the given breakends
	 * @param b1
	 * @param b2
	 * @return
	 */
	public static BreakendSummary overlapOf(BreakendSummary b1,
			BreakendSummary b2) {
		if (b1.referenceIndex != b2.referenceIndex) return null;
		if (b1.direction != b2.direction) return null;
		int start = Math.max(b1.start, b2.start);
		int end = Math.min(b1.end, b2.end);
		if (start > end) return null;
		BreakendSummary result = new BreakendSummary(b1.referenceIndex, b1.direction, start, end);
		return result;
	}
	/**
	 * Extends the breakend location interval
	 * @param expandBy number of bases to extend bounds by
	 * @param dictionary sequence dictionary
	 * @return breakend with bounds expanded on both sides
	 */
	public BreakendSummary expandBounds(int expandBy, SAMSequenceDictionary dictionary) {
		return new BreakendSummary(referenceIndex, direction, Math.max(1, start - expandBy), Math.min(dictionary.getSequence(referenceIndex).getSequenceLength(), end + expandBy));
	}
	/**
	 * Reduces size of the breakend location interval
	 * @param by bases to reduce each side of the bounds by 
	 * @param dictionary sequence dictionary
	 * @return breakend with bounds reduced
	 */
	public BreakendSummary compressBounds(int by) {
		return compressBounds(by, RoundingMode.DOWN);
	}
	protected BreakendSummary compressBounds(int by, RoundingMode roundingMode) {
		if (end - start + 1 <= 2 * by) {
			int centre = IntMath.divide(end + start, 2, roundingMode);
			return new BreakendSummary(referenceIndex, direction, centre, centre);
		}
		return new BreakendSummary(referenceIndex, direction, start + by, end - by);
	}
	protected static String toString(int referenceIndex, int start, int end, ProcessingContext processContext) {
		if (processContext != null && referenceIndex >= 0 && referenceIndex < processContext.getDictionary().size()) {
			return String.format("%s:%d-%d", processContext.getDictionary().getSequence(referenceIndex).getSequenceName(), start, end);
		}
		return String.format("(%d):%d-%d", referenceIndex, start, end);
	}
	protected static String toString(BreakendDirection direction, int referenceIndex, int start, int end, ProcessingContext processContext) {
		if (direction == BreakendDirection.Forward) return toString(referenceIndex, start, end, processContext) + ">";
		return "<" + toString(referenceIndex, start, end, processContext);
	}
	/**
	 * Determines whether this given breakend is valid for the given reference
	 * @param dictionary
	 * @return
	 */
	public boolean isValid(SAMSequenceDictionary dictionary) {
		return referenceIndex >= 0 && referenceIndex < dictionary.size()
				&& start <= end
				&& start > 0 && end <= dictionary.getSequence(referenceIndex).getSequenceLength();
	}
	@Override
	public String toString() {
		return toString(direction, referenceIndex, start, end, null);
	}
	public String toString(ProcessingContext processContext) {
		return toString(direction, referenceIndex, start, end, processContext);
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((direction == null) ? 0 : direction.hashCode());
		result = prime * result + end;
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
}
