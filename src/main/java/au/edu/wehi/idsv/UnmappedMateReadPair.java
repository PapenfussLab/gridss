package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public class UnmappedMateReadPair extends NonReferenceReadPair {
	protected UnmappedMateReadPair(SAMRecord local, SAMRecord remote, SAMEvidenceSource source, ProcessingContext processContext) {
		super(local, remote, source, processContext);
	}
	@Override
	public String toString() {
		return String.format("UM %s MQ=%d RN=%s", getBreakendSummary(), getLocalMapq(), getLocalledMappedRead().getReadName());
	}
}
