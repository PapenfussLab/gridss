package au.edu.wehi.idsv;

import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MathUtil;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.math.IntMath;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;

import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;

/**
 * Positional location of a breakpoint that is consistent with the given evidence
 * @author Daniel Cameron
 *
 */
public class BreakendSummary {
	/**
	 * Nominal position adjacent to breakpoint in 1-based genomic coordinate
	 */
	public final int nominal;
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
	public BreakendSummary(int referenceIndex, BreakendDirection direction, int nominal) {
		this(referenceIndex, direction, nominal, nominal, nominal);
	}
	public BreakendSummary(int referenceIndex, BreakendDirection direction, int nominal, int start, int end) {
		if (referenceIndex < 0) {
			throw new IllegalArgumentException("Reference index must be valid");
		}
		if (end < start) {
			throw new IllegalArgumentException("end must be at or after start");
		}
		if (nominal < start || nominal > end) {
			throw new IllegalArgumentException("nominal must be within [start, end] interval");
		}
		this.referenceIndex = referenceIndex;
		this.direction = direction;
		this.start = start;
		this.nominal = nominal;
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
				IntervalUtil.overlapsClosed(this.start, this.end, loc.start, loc.end);
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
		int nominal = boundedNominal(MathUtil.average(b1.nominal, b2.nominal), start, end); 
		BreakendSummary result = new BreakendSummary(b1.referenceIndex, b1.direction, nominal, start, end);
		return result;
	}
	/**
	 * Returns the nominal position whilst ensuring it is within the given bounds
	 * @param start
	 * @param nominal
	 * @param end
	 * @return
	 */
	private static int boundedNominal(int nominal, int start, int end) {
		return Math.min(Math.max(nominal, start), end);
	}
	/**
	 * Extends the breakend location interval
	 * @param expandBy number of bases to extend bounds by
	 * @param dictionary sequence dictionary
	 * @return breakend with bounds expanded on both sides
	 */
	public BreakendSummary expandBounds(int expandBy) {
		return adjustBounds(-expandBy, expandBy);
	}
	/**
	 * Extends the breakend location interval
	 * @param expandBy number of bases to extend bounds by
	 * @param dictionary sequence dictionary
	 * @return breakend with bounds expanded on both sides
	 */
	private BreakendSummary adjustBounds(int startAdjustment, int endAdjustment) {
		return adjustBounds(startAdjustment, endAdjustment, RoundingMode.DOWN);
	}
	private BreakendSummary adjustBounds(int startAdjustment, int endAdjustment, RoundingMode roundingMode) {
		int newStart = start + startAdjustment;
		int newEnd = end + endAdjustment;
		if (newStart > newEnd) {
			newStart = IntMath.divide(newStart + newEnd, 2, roundingMode);
			newEnd = newStart;
		}
		int newNominal = boundedNominal(nominal, newStart, newEnd);
		return new BreakendSummary(referenceIndex, direction, newNominal, newStart, newEnd);
	}
	/**
	 * Adjusts the position of the breakend anchor
	 * @param intoBreakendAnchor number of anchored bases to expand the breakend bounds into
	 * @param awayFromBreakendAnchor number of additional anchored bases to consider
	 * @return
	 */
	public BreakendSummary adjustPosition(int intoBreakendAnchor, int awayFromBreakendAnchor, boolean adjustNominal) {
		if (adjustNominal && (intoBreakendAnchor < 0 || awayFromBreakendAnchor < 0)) {
			throw new IllegalArgumentException("Cannot make adjustment that causes the nominal position to no longer be valid");
		}
		if (intoBreakendAnchor == 0 && awayFromBreakendAnchor == 0) return this;
		int anchorShift = (-intoBreakendAnchor + awayFromBreakendAnchor) / 2;
		int startAdjust = -intoBreakendAnchor;
		int endAdjust = awayFromBreakendAnchor;
		if (direction == BreakendDirection.Backward) {
			startAdjust = -awayFromBreakendAnchor;
			endAdjust = intoBreakendAnchor;
			anchorShift *= -1;
		}
		if (!adjustNominal) {
			anchorShift = 0;
		}
		return new BreakendSummary(referenceIndex, direction, nominal + anchorShift, start + startAdjust, end + endAdjust);
	}
	/**
	 * Reduces size of the breakend location interval
	 * @param by bases to reduce each side of the bounds by 
	 * @param dictionary sequence dictionary
	 * @return breakend with bounds reduced
	 */
	public BreakendSummary compressBounds(int by) {
		return adjustBounds(by, -by);
	}
	protected BreakendSummary compressBounds(int by, RoundingMode roundingMode) {
		return adjustBounds(by, -by, roundingMode);
	}
	/**
	 * Gets the nominal position of the variant for variant calling purposes
	 * @return position to call
	 */
	public BreakendSummary getNominalPosition() {
		return new BreakendSummary(referenceIndex, direction, nominal, nominal, nominal);
	}
	protected static String toString(int referenceIndex, int nominal, int start, int end, GenomicProcessingContext processContext) {
		StringBuilder sb = new StringBuilder();
		if (processContext != null && referenceIndex >= 0 && referenceIndex < processContext.getDictionary().size()) {
			sb.append(processContext.getDictionary().getSequence(referenceIndex).getSequenceName());
		} else {
			sb.append('(');
			sb.append(referenceIndex);
			sb.append(')');
		}
		sb.append(':');
		sb.append(nominal);
		if (start != nominal || end != nominal) {
			sb.append('(');
			sb.append(start);
			sb.append('-');
			sb.append(end);
			sb.append(')');
		}
		return sb.toString();
	}
	protected static String toString(BreakendDirection direction, int referenceIndex, int nominal, int start, int end, GenomicProcessingContext processContext) {
		if (direction == BreakendDirection.Forward) return toString(referenceIndex, nominal, start, end, processContext) + ">";
		return "<" + toString(referenceIndex, nominal, start, end, processContext);
	}
	/**
	 * Determines whether this given breakend is valid for the given reference
	 * @param dictionary
	 * @return
	 */
	public boolean isValid(SAMSequenceDictionary dictionary) {
		if (dictionary == null) {
			throw new NullPointerException("Missing required .dict sequence dictionary");
		}
		return referenceIndex >= 0 && referenceIndex < dictionary.size()
				&& start <= end
				&& start > 0 && end <= dictionary.getSequence(referenceIndex).getSequenceLength();
	}
	/**
	 * Adjusts the breakend positions such that the call variant positions are valid position
	 * for the contig in question.
	 * @param dictionary
	 * @return adjusted positions
	 */
	public BreakendSummary asValidFor(SAMSequenceDictionary dictionary) {
		if (isValid(dictionary)) return this;
		int newStart = within(start, 1, dictionary.getSequence(referenceIndex).getSequenceLength());
		int newEnd = within(end, 1, dictionary.getSequence(referenceIndex).getSequenceLength());
		return new BreakendSummary(referenceIndex, direction, boundedNominal(nominal, newStart, newEnd), newStart, newEnd);
	}
	private static int within(int x, int min, int max) {
		return Math.min(Math.max(x, min), max);
	}
	/**
	 * Gets a minimal CIGAR representing this breakend.
	 * @return
	 */
	public List<CigarElement> getCigarRepresentation() {
		List<CigarElement> breakcigar = new ArrayList<CigarElement>(3);
		if (start == end) {
			breakcigar.add(new CigarElement(1, CigarOperator.X));
		} else if (end - start == 1) {
			breakcigar.add(new CigarElement(2, CigarOperator.X));
		} else {
			breakcigar.add(new CigarElement(1, CigarOperator.X));
			breakcigar.add(new CigarElement(end - start - 1, CigarOperator.N));
			breakcigar.add(new CigarElement(1, CigarOperator.X));
		}
		return breakcigar;
	}
	@Override
	public String toString() {
		return toString(direction, referenceIndex, nominal, start, end, null);
	}
	public String toString(GenomicProcessingContext processContext) {
		return toString(direction, referenceIndex, nominal, start, end, processContext);
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((direction == null) ? 0 : direction.hashCode());
		result = prime * result + nominal;
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
		if (nominal != other.nominal)
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
			        .compare(o1.nominal, o2.nominal)
			        .result();
		  }
	};
	public static Ordering<BreakendSummary> ByEndStart = new Ordering<BreakendSummary>() {
		public int compare(BreakendSummary o1, BreakendSummary o2) {
			  return ComparisonChain.start()
			        .compare(o1.referenceIndex, o2.referenceIndex)
			        .compare(o1.end, o2.end)
			        .compare(o1.start, o2.start)
			        .compare(o1.nominal, o2.nominal)
			        .result();
		  }
	};
}
