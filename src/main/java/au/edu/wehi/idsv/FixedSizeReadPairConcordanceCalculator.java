package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;

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
	public int minConcordantFragmentSize() {
		return minFragmentSize;
	}
	@Override
	public boolean isConcordant(SAMRecord read1, SAMRecord read2) {
		return super.isConcordant(read1, read2)
				&& isConcordantFragmentSize(read1, read2);
	}
	@Override
	public boolean isConcordant(SAMRecord read) {
		return isConcordant(read, null);
	}
	private boolean isConcordantFragmentSize(SAMRecord read1, SAMRecord read2) {
		int fragSize;
		if (read2 == null) {
			fragSize = SAMRecordUtil.estimateFragmentSize(read1, PairOrientation.FR);
		} else {
			fragSize = SAMRecordUtil.calculateFragmentSize(read1, read2, PairOrientation.FR);
		}
		return minFragmentSize <= fragSize && fragSize <= maxFragmentSize;
	}
}
