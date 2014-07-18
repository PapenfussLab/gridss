package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public class UnmappedMateReadPair extends NonReferenceReadPair {
	@Override
	public float getPhredLogLikelihoodRatio() {
		// TODO: real metrics
		return getLocalMapq();
	}
	protected UnmappedMateReadPair(SAMRecord local, SAMRecord remote, SAMEvidenceSource source) {
		super(local, remote, source);
	}
}
