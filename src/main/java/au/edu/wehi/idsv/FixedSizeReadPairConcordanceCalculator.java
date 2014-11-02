package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class FixedSizeReadPairConcordanceCalculator extends ReadPairConcordanceCalculator {
	private final int minFragmentSize;
	private final int maxFragmentSize;
	public FixedSizeReadPairConcordanceCalculator(int minFragmentSize, int maxFragmentSize) {
		this.minFragmentSize = Math.max(1, minFragmentSize);
		this.maxFragmentSize = maxFragmentSize;
	}
	@Override
	public int maxConcordantFragmentSize() {
		return maxFragmentSize;
	}
	@Override
	public boolean isConcordant(SAMRecord local) {
		return local.getReadPairedFlag()
				&& !local.getReadUnmappedFlag()
				&& !local.getMateUnmappedFlag()
				&& local.getReferenceIndex() == local.getMateReferenceIndex()
				// (assumes FR)
				&& local.getReadNegativeStrandFlag() != local.getMateNegativeStrandFlag()
				&& isConcordantFragmentSize(local);
	}
	private boolean isConcordantFragmentSize(SAMRecord local) {
		int fragSize = SAMRecordUtil.estimateFragmentSize(local);
		return minFragmentSize <= fragSize && fragSize <= maxFragmentSize;
	}
}
