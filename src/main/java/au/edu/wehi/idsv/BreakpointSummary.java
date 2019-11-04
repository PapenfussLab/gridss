package au.edu.wehi.idsv;

import au.edu.wehi.idsv.debruijn.positional.SupportNodeIterator;
import au.edu.wehi.idsv.util.IntervalUtil;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;

import java.math.RoundingMode;

/**
 * Positional locations on source and target chromosomes
 * of a breakpoint that is consistent with the given evidence
 * @author Daniel Cameron
 *
 */
public class BreakpointSummary extends BreakendSummary {
	private static final Log log = Log.getInstance(BreakpointSummary.class);
	/**
	 * Nominal position of breakpoint is immediately after this 1-based genomic coordinate on the destination contig
	 */
	public final int nominal2;
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
	public BreakpointSummary(int referenceIndex1, BreakendDirection direction1, int nominal1,
			int referenceIndex2, BreakendDirection direction2, int nominal2) {
		this(referenceIndex1, direction1, nominal1, nominal1, nominal1,
				referenceIndex2, direction2, nominal2, nominal2, nominal2);
	}
	public BreakpointSummary(int referenceIndex1, BreakendDirection direction1, int nominal1, int start1, int end1,
			int referenceIndex2, BreakendDirection direction2, int nominal2, int start2, int end2) {
		super(referenceIndex1, direction1, nominal1, start1, end1);
		if (referenceIndex2 < 0) {
			throw new IllegalArgumentException("Reference index must be valid");
		}
		this.nominal2 = nominal2;
		this.start2 = start2;
		this.end2 = end2;
		this.referenceIndex2 = referenceIndex2;
		this.direction2 = direction2;
		if (end2 < start2) {
			throw new IllegalArgumentException("end must be at or after start");
		}
		if (nominal2 < start2 || nominal2 > end2) {
			throw new IllegalArgumentException("nominal must be within [start, end] interval");
		}
	}
	public BreakpointSummary(BreakendSummary local, BreakendSummary remote) {
		this(local.referenceIndex, local.direction, local.nominal, local.start, local.end, 
				remote.referenceIndex, remote.direction, remote.nominal, remote.start, remote.end);
	}
	public BreakendSummary localBreakend() {
		return new BreakendSummary(referenceIndex, direction, nominal, start, end);
	}
	public BreakendSummary remoteBreakend() {
		return new BreakendSummary(referenceIndex2, direction2, nominal2, start2, end2);
	}
	/**
	 * Determines whether this breakend is the higher of the two breakends
	 * @return true if this breakend starts after the remote breakend
	 */
	public boolean isHighBreakend() {
		if (referenceIndex != referenceIndex2) return referenceIndex > referenceIndex2;
		if (start != start2) return start > start2;
		if (end != end2) return end > end2;
		return nominal > nominal2;
	}
	public BreakendSummary highBreakend() {
		if (isHighBreakend()) return localBreakend();
		return remoteBreakend();
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
	public BreakendSummary lowBreakend() {
		if (isLowBreakend()) return localBreakend();
		return remoteBreakend();
	}
	/**
	 * Returns the other end of this breakpoint
	 * @return breakpoint with ends swapped
	 */
	public BreakpointSummary remoteBreakpoint() {
		return new BreakpointSummary(referenceIndex2, direction2, nominal2, start2, end2, referenceIndex, direction, nominal, start, end);
	}
	@Override
	public String toString() {
		return String.format("%s %s", toString(direction, referenceIndex, nominal, start, end, null), toString(direction2, referenceIndex2, nominal2, start2, end2, null));
	}
	@Override
	public String toString(GenomicProcessingContext processContext) {
		return String.format("%s %s", toString(direction, referenceIndex, nominal, start, end, processContext), toString(direction2, referenceIndex2, nominal2, start2, end2, processContext));
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
	/**
	 * Adjusts the position of the breakpoint anchor relative to the local breakend position.
	 * Adjustments are made to both breakends to conserve the event size.
	 * @param intoBreakendAnchor number of additional bases to allow to be allocated to the remote breakend 
	 * @param awayFromBreakendAnchor number of additional bases to allow to be allocated to the local breakend
	 * @return expanded breakpoint
	 */
	@Override
	public BreakpointSummary adjustPosition(int intoBreakendAnchor, int awayFromBreakendAnchor, boolean adjustNominal) {
		if (intoBreakendAnchor == 0 && awayFromBreakendAnchor == 0) return this;
		return new BreakpointSummary(super.adjustPosition(intoBreakendAnchor, awayFromBreakendAnchor, adjustNominal),
				remoteBreakend().adjustPosition(awayFromBreakendAnchor, intoBreakendAnchor, adjustNominal));
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
	public boolean couldBeReferenceAllele(boolean isPrecise) {
		if (!isPrecise) {
			return couldBeDeletionOfSize(0, 0);
		} else {
			return new BreakpointSummary(referenceIndex, direction, nominal, referenceIndex2, direction2, nominal2)
				.couldBeReferenceAllele(false);
		}
	}
	/**
	 * Gets the nominal position of the variant for variant calling purposes
	 * @return position to call
	 */
	@Override
	public BreakpointSummary getNominalPosition() {
		return new BreakpointSummary(
				referenceIndex, direction, nominal, nominal, nominal,
				referenceIndex2, direction2, nominal2, nominal2, nominal2);
	}
	/**
	 * Adjusts the breakend positions such that the call variant positions are valid position
	 * for the contigs in question.
	 * @param dictionary
	 * @return adjusted positions
	 */
	public BreakpointSummary asValidFor(SAMSequenceDictionary dictionary) {
		if (isValid(dictionary)) return this;
		return new BreakpointSummary(super.asValidFor(dictionary), remoteBreakend().asValidFor(dictionary));
	}
	/**
	 * Determines whether this given breakpoint is valid for the given reference
	 * @param dictionary
	 * @return
	 */
	public boolean isValid(SAMSequenceDictionary dictionary) {
		return super.isValid(dictionary) &&
				referenceIndex2 >= 0 && referenceIndex2 < dictionary.size()
				&& start2 <= end2
				&& start2 > 0 && end2 <= dictionary.getSequence(referenceIndex2).getSequenceLength();
	}
	/**
	 * Determines the size of the simplest event (deletion, inversion, tandem duplication)
	 * this breakpoint contributes to.
	 * @return event size, null for interchromosomal breakpoints
	 */
	public Integer getEventSize() {
		if (referenceIndex != referenceIndex2) return null;
		if (direction == direction2) {
			return Math.abs(start2 - start);
		} else if ((start < start2 && direction == BreakendDirection.Forward) || 
				(start2 < start && direction == BreakendDirection.Backward)) {
			return Math.abs(start2 - start) - 1;
		} else {
			return Math.abs(start2 - start) + 1;
		}
	}
	@Override
	public BreakpointSummary centreAligned() {
		int addToLocal;
		int addToRemote;
		double targetPosition = ((double)start + end) / 2;
		if (direction != direction2) {
			addToLocal = (int)Math.floor(targetPosition) - nominal;
			addToRemote = addToLocal;
		} else {
			if (isLowBreakend()) {
				addToLocal = (int)Math.floor(targetPosition) - nominal;
			} else {
				addToLocal = (int)Math.ceil(targetPosition) - nominal;
			}
			addToRemote = -addToLocal;
		}
		try {
			return new BreakpointSummary(
					referenceIndex, direction, nominal + addToLocal, start, end,
					referenceIndex2, direction2, nominal2 + addToRemote, start2, end2);
		} catch (IllegalArgumentException e) {
			log.debug("Adjusted %s out of interval bounds.", this.toString());
			return this;
		}
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
			        .compare(o1.nominal, o2.nominal)
			        .compare(o1.nominal2, o2.nominal2)
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