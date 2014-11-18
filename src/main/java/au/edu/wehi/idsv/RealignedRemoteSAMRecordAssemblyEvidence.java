package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public class RealignedRemoteSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence implements RemoteEvidence {
	public RealignedRemoteSAMRecordAssemblyEvidence(ProcessingContext processContext, AssemblyEvidenceSource source, SAMRecord assembly, SAMRecord realigned) {
		super(processContext, source, assembly, realigned);
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return getBreakendSummary().remoteBreakpoint();
	}
}
