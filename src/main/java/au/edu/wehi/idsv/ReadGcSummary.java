package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.reference.ReferenceSequence;

public class ReadGcSummary {
	/**
	 * GC content assumed when GC cannot be calculated (e.g. reference bases are all N).
	 */
	public static final double UNDEFINED_GC = 50;
	public final int referenceIndex;
	public final int fragmentStart;
	public final int fragmentEnd;
	public final double gcPercentage;
	public ReadGcSummary(SAMRecord record, ReferenceSequence refSeq, int defaultFragmentSize, ReadPairConcordanceCalculator rpcc) {
		if (record.getReadUnmappedFlag()) {
			throw new IllegalArgumentException("Read must be mapped.");
		}
		referenceIndex = record.getReferenceIndex();
		int fragmentSize = getFragmentSize(record, defaultFragmentSize, rpcc);
    	if (record.getReadNegativeStrandFlag()) {
    		this.fragmentEnd = record.getAlignmentEnd();
    		this.fragmentStart = fragmentEnd - fragmentSize + 1;
    	} else {
    		this.fragmentStart = record.getAlignmentStart();
    		this.fragmentEnd = fragmentStart + fragmentSize - 1;
    	}
    	this.gcPercentage = getReferenceGCPercentage(fragmentStart -1, fragmentEnd, refSeq);    	
	}
	private static double getReferenceGCPercentage(int zeroBasedStartInclusive, int zeroBasedEndExclusive, ReferenceSequence refSeq) {
    	byte[] ref = refSeq.getBases();
    	int gcCount = 0;
    	int atCount = 0;
    	for (int i = zeroBasedStartInclusive; i < zeroBasedEndExclusive; i++) {
    		if (ref[i] == 'G' || ref[i] == 'C' || ref[i] == 'g' || ref[i] == 'c') {
    			gcCount++;
    		}
    		if (ref[i] == 'A' || ref[i] == 'T' || ref[i] == 'a' || ref[i] == 't') {
    			atCount++;
    		}
    	}
    	if (gcCount + atCount == 0) {
    		return UNDEFINED_GC;
    	}
    	return gcCount / (double)(gcCount + atCount);
    }
	private static int getFragmentSize(SAMRecord samRecord, int defaultFragmentSize, ReadPairConcordanceCalculator rpcc) {
    	int fragmentSize = Math.max(defaultFragmentSize, samRecord.getReadLength());
    	if (rpcc.isConcordant(samRecord) && samRecord.getMateAlignmentStart() != 0) {
    		fragmentSize = SAMRecordUtil.estimateFragmentSize(samRecord, PairOrientation.FR); 
    	}
    	return fragmentSize;
    }
}
