package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

/**
 * Soft clip split read in which the realigned soft clip maps to a local coordinate
 * 
 * @author Daniel Cameron
 *
 */
public class RealignedRemoteSoftClipEvidence extends RealignedSoftClipEvidence {
	private final BreakpointSummary location;
	protected RealignedRemoteSoftClipEvidence(ProcessingContext processContext,
			SAMEvidenceSource source, BreakendDirection direction,
			SAMRecord record, SAMRecord realigned)
			throws CloneNotSupportedException {
		super(processContext, source, direction, record, realigned);
		this.location = super.getBreakendSummary().remoteBreakpoint();
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return location;
	}
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return "Remote" + super.toString();
	}
	// TODO: should we flip getLocal*() and getRemote*() so local = here
	// or leave unchanged so local = read anchor (which in this case is at the remote breakend)
	// or special case all handling?
}
