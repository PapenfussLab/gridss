package au.edu.wehi.socrates;

import htsjdk.samtools.SAMRecord;

public class RealignedBreakpoint {
	private final SAMRecord realigned;
	private final BreakpointSummary summary;
	public RealignedBreakpoint(BreakendSummary local, SAMRecord realigned) {
		if (realigned.getReadUnmappedFlag()) throw new IllegalArgumentException("Realigned read is not mapped");
		this.realigned = realigned;
		// TODO: refactor SC and assembly realignment logic into here
	}
	public int getInsertedSequenceLength() {
		return getInsertedSequence().length();
	}
	public BreakpointSummary getBreakpointSummary() {
		return summary;
	}
	/**
	 * Returns the untemplated breakpoint sequence not mapped to either the local
	 * or realigned mapping location.
	 * @return
	 */
	public String getInsertedSequence() {
	}
	public SAMRecord getSAMRecord() {
		return realigned;
	}
}
