package au.edu.wehi.socrates;

/**
 * Positional locations on source and target chromosomes
 * of a breakpoint that is consistent with the given evidence
 * @author Daniel Cameron
 *
 */
public class BreakpointInterval extends BreakpointLocation {
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
	public final BreakpointDirection direction2;
	public BreakpointInterval(int referenceIndex1, BreakpointDirection direction1, int start1, int end1,
			int referenceIndex2, BreakpointDirection direction2, int start2, int end2, double qual) {
		super(referenceIndex1, direction1, start1, end1, qual);
		this.start2 = start2;
		this.end2 = end2;
		this.referenceIndex2 = referenceIndex2;
		this.direction2 = direction2;
	}
	public BreakpointInterval(BreakpointLocation local, BreakpointLocation remote, double qual) {
		this(local.referenceIndex, local.direction, local.start, local.end,
				remote.referenceIndex, remote.direction, remote.start, remote.end,
				qual);
	}
	public BreakpointInterval(BreakpointInterval interval, double qual) {
		this(interval.referenceIndex, interval.direction, interval.start, interval.end,
				interval.referenceIndex2, interval.direction2, interval.start2, interval.end2,
				qual);
	}
	@Override
	public String toString() {
		return String.format("%s %s %.1f", toString(direction, referenceIndex, start, end), toString(direction2, referenceIndex2, start2, end2), qual);
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
		BreakpointInterval other = (BreakpointInterval) obj;
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
}