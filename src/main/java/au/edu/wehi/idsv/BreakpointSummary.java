package au.edu.wehi.idsv;

import java.math.RoundingMode;

import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.math.IntMath;

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
	 * Determines whether this breakend is the higher of the two breakends
	 * @return true if this breakend starts after the remote breakend
	 */
	public boolean isHighBreakend() {
		if (referenceIndex != referenceIndex2) return referenceIndex > referenceIndex2;
		if (start != start2) return start > start2;
		return end > end2;
	}
	/**
	 * Determines whether this breakend is the lower of the two breakends
	 * @return true if this breakend starts before the remote breakend
	 */
	public boolean isLowBreakend() {
		if (referenceIndex != referenceIndex2) return referenceIndex < referenceIndex2;
		if (start != start2) return start < start2;
		return end < end2;
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
	public String toString(GenomicProcessingContext processContext) {
		return String.format("%s %s", toString(direction, referenceIndex, start, end, processContext), toString(direction2, referenceIndex2, start2, end2, processContext));
	}
	/**
	 * Restricts overlap to require both local and remote breakends to overlap
	 */
	@Override
	public boolean overlaps(BreakendSummary loc) {
		if (loc instanceof BreakpointSummary) return overlaps((BreakpointSummary)loc);
		return super.overlaps(loc);
	}
	public boolean overlaps(BreakpointSummary loc) {
		return breakendOverlaps(loc) && remoteBreakendOverlaps((BreakpointSummary)loc);
	}
	protected boolean remoteBreakendOverlaps(BreakpointSummary loc) {
		return this.referenceIndex2 == loc.referenceIndex2 &&
				this.direction2 == loc.direction2 &&
				IntervalUtil.overlapsClosed(this.start2, this.end2, loc.start2, loc.end2);
	}
	@Override
	public BreakpointSummary expandBounds(int expandBy) {
		return new BreakpointSummary(localBreakend().expandBounds(expandBy), remoteBreakend().expandBounds(expandBy));
	}
	@Override
	public BreakpointSummary compressBounds(int by) {
		BreakendSummary local = localBreakend();
		BreakendSummary remote = remoteBreakend();
		int order = BreakendSummary.ByStartEnd.compare(local, remote);
		if (order < 0) {
			return new BreakpointSummary(local.compressBounds(by, RoundingMode.DOWN), remote.compressBounds(by, RoundingMode.UP));
		} else  if (order == 0) {
			return new BreakpointSummary(local.compressBounds(by, RoundingMode.DOWN), remote.compressBounds(by, RoundingMode.DOWN));
		} else {
			return new BreakpointSummary(local.compressBounds(by, RoundingMode.UP), remote.compressBounds(by, RoundingMode.DOWN));
		}
	}
	public boolean couldBeDeletionOfSize(int minSize, int maxSize) {
		if (referenceIndex != referenceIndex2 || direction == direction2) return false;
		int fwdStart, fwdEnd, bwdStart, bwdEnd;
		if (direction == BreakendDirection.Forward) {
			fwdStart = start;
			fwdEnd = end;
			bwdStart = start2;
			bwdEnd = end2;
		} else {
			bwdStart = start;
			bwdEnd = end;
			fwdStart = start2;
			fwdEnd = end2;
		}
		int bpMinSize = bwdStart - fwdEnd - 1;
		int bpMaxSize = bwdEnd - fwdStart - 1;
		return IntervalUtil.overlapsClosed(bpMinSize, bpMaxSize, minSize, maxSize);
	}
	/**
	 * Determines whether this breakpoint could just be the reference allele
	 * @return
	 */
	public boolean couldBeReferenceAllele() {
		return couldBeDeletionOfSize(0, 0);
	}
	/**
	 * Gets the nominal position of the variant for variant calling purposes
	 * @return position to call
	 */
	public BreakpointSummary getCallPosition() {
		int callPos;
		int remoteCallPos;
		// round the called position of the lower breakend down
		if (BreakendSummary.ByStartEnd.compare(localBreakend(), remoteBreakend()) > 0) {
			callPos = IntMath.divide(start + end, 2, RoundingMode.CEILING);
			remoteCallPos = IntMath.divide(start2 + end2, 2, RoundingMode.FLOOR);
		} else {
			callPos = IntMath.divide(start + end, 2, RoundingMode.FLOOR);
			remoteCallPos = IntMath.divide(start2 +end2, 2, RoundingMode.CEILING);
		}
		return new BreakpointSummary(referenceIndex, direction, callPos, callPos, referenceIndex2, direction2, remoteCallPos, remoteCallPos);
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
	/**
	 * Orders breakpoints in ascending order of the genomic location of both breakends 
	 */
	public static Ordering<BreakpointSummary> ByLowHigh = new Ordering<BreakpointSummary>() {
		public int compare(BreakpointSummary o1, BreakpointSummary o2) {
			return ByStartStart2EndEnd2.compare(o1.isHighBreakend() ? o1.remoteBreakpoint() : o1, o2.isHighBreakend() ? o2.remoteBreakpoint() : o2);
		  }
	};
}