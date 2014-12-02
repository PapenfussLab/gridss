package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

/**
 * Positional locations on source and target chromosomes
 * of a breakpoint that is consistent with the given evidence
 * @author Daniel Cameron
 *
 */
public class BreakpointSummary extends BreakendSummary {
	/**
	 * First possible position of breakpoint is immediately after this 1-based genomic coordinate on the destination contig
	 */
	public final int start2;
	/**
	 * Last possible position of breakpoint is immediately after this 1-based genomic coordinate on the destination contig
	 */
	public final int end2;
	/**
	 * Destination contig of breakpoint
	 */
	public final int referenceIndex2;
	/**
	 * Breakpoint is in the given direction on the destination contig
	 */
	public final BreakendDirection direction2;
	public BreakpointSummary(int referenceIndex1, BreakendDirection direction1, int start1, int end1,
			int referenceIndex2, BreakendDirection direction2, int start2, int end2) {
		super(referenceIndex1, direction1, start1, end1);
		if (referenceIndex2 < 0) {
			throw new IllegalArgumentException("Reference index must be valid");
		}
		this.start2 = start2;
		this.end2 = end2;
		this.referenceIndex2 = referenceIndex2;
		this.direction2 = direction2;
	}
	public BreakpointSummary(BreakendSummary local, BreakendSummary remote) {
		this(local.referenceIndex, local.direction, local.start, local.end,
				remote.referenceIndex, remote.direction, remote.start, remote.end);
	}
	public BreakendSummary localBreakend() {
		return new BreakendSummary(referenceIndex, direction, start, end);
	}
	public BreakendSummary remoteBreakend() {
		return new BreakendSummary(referenceIndex2, direction2, start2, end2);
	}
	/**
	 * Returns the other end of this breakpoint
	 * @return breakpoint with ends swapped
	 */
	public BreakpointSummary remoteBreakpoint() {
		return new BreakpointSummary(referenceIndex2, direction2, start2, end2, referenceIndex, direction, start, end);
	}
	@Override
	public String toString() {
		return String.format("%s %s", toString(direction, referenceIndex, start, end, null), toString(direction2, referenceIndex2, start2, end2, null));
	}
	@Override
	public String toString(ProcessingContext processContext) {
		return String.format("%s %s", toString(direction, referenceIndex, start, end, processContext), toString(direction2, referenceIndex2, start2, end2, processContext));
	}
	/**
	 * Restricts overlap to require both local and remote breakends to overlap
	 */
	@Override
	public boolean overlaps(BreakendSummary loc) {
		return super.overlaps(loc) &&
				(!(loc instanceof BreakpointSummary) ||
						remoteBreakend().breakendOverlaps(((BreakpointSummary)loc).remoteBreakend()));
	}
	@Override
	public BreakpointSummary expandBounds(int expandBy, SAMSequenceDictionary dictionary) {
		return new BreakpointSummary(localBreakend().expandBounds(expandBy, dictionary), remoteBreakend().expandBounds(expandBy, dictionary));
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result
				+ ((direction2 == null) ? 0 : direction2.hashCode());
		result = prime * result + end2;
		result = prime * result + referenceIndex2;
		result = prime * result + start2;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		BreakpointSummary other = (BreakpointSummary) obj;
		if (direction2 != other.direction2)
			return false;
		if (end2 != other.end2)
			return false;
		if (referenceIndex2 != other.referenceIndex2)
			return false;
		if (start2 != other.start2)
			return false;
		return true;
	}
	public static Ordering<BreakpointSummary> ByStartStart2EndEnd2 = new Ordering<BreakpointSummary>() {
		public int compare(BreakpointSummary o1, BreakpointSummary o2) {
			  return ComparisonChain.start()
			        .compare(o1.referenceIndex, o2.referenceIndex)
			        .compare(o1.start, o2.start)
			        .compare(o1.referenceIndex2, o2.referenceIndex2)
			        .compare(o1.start2, o2.start2)
			        .compare(o1.end, o2.end)
			        .compare(o1.end2, o2.end2)
			        .result();
		  }
	};
}