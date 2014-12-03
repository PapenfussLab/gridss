package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public class RealignedRemoteSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence implements RemoteEvidence {
	private final BreakpointSummary bs;
	public RealignedRemoteSAMRecordAssemblyEvidence(ProcessingContext processContext, AssemblyEvidenceSource source, SAMRecord assembly, SAMRecord realigned) {
		super(processContext, source, assembly, realigned);
		this.bs = super.getBreakendSummary().remoteBreakpoint();
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return bs;
	}
	@Override
	public String toString() {
		return "R" + super.toString();
	}
}
