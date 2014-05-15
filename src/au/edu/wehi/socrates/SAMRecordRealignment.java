package au.edu.wehi.socrates;

import htsjdk.samtools.SAMRecord;

/**
 * Helper class for calculating breakpoint properties from a breakend and realignment
 * SAMRecord of the breakpoint sequence 
 * @author Daniel Cameron
 *
 */
public class SAMRecordRealignment {
	private final BreakendSummary breakend;
	private final SAMRecord realignment;
	public SAMRecordRealignment(BreakendSummary breakend, SAMRecord realignment) {
		this.breakend = breakend;
		this.realignment = realignment;
	}
	public BreakpointSummary getBreakpointSummary() {
	}
	public int getUntemplatedSequenceLength() {
	}
}
