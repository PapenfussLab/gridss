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
	public final long start2;
	/**
	 * Last possible position of breakpoint is immediately after this 1-based genomic coordinate on the destination contig
	 */
	public final long end2;
	/**
	 * Destination contig of breakpoint
	 */
	public final int referenceIndex2;
	/**
	 * Breakpoint is in the given direction on the destination contig
	 */
	public final BreakpointDirection direction2;
	public BreakpointInterval(int referenceIndex1, BreakpointDirection direction1, long start1, long end1,
			int referenceIndex2, BreakpointDirection direction2, long start2, long end2, double qual) {
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
}