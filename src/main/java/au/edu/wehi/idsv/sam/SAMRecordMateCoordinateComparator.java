package au.edu.wehi.idsv.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;

/**
 * Comparator for sorting SAMRecords by mate coordinate.
 *
 */
public class SAMRecordMateCoordinateComparator implements SAMRecordComparator {
	public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
		int cmp = fileOrderCompare(samRecord1, samRecord2);
		// Note: secondary sort order does not match SAMRecordCoordinateComparator
		if (cmp != 0) return cmp;
		cmp = samRecord1.getReadName().compareTo(samRecord2.getReadName());
		if (cmp != 0) return cmp;
		cmp = compareInts(samRecord1.getFlags(), samRecord2.getFlags());
		if (cmp != 0) return cmp;
		cmp = compareInts(samRecord1.getMappingQuality(), samRecord2.getMappingQuality());
		if (cmp != 0) return cmp;
		cmp = compareInts(samRecord1.getInferredInsertSize(), samRecord2.getInferredInsertSize());
		return cmp;
    }

    private int compareInts(int i1, int i2) {
        if (i1 < i2) return -1;
        else if (i1 > i2) return 1;
        else return 0;
    }

    /**
     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
     */
	public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        final int refIndex1 = samRecord1.getMateReferenceIndex();
        final int refIndex2 = samRecord2.getMateReferenceIndex();
        if (refIndex1 == -1) {
            return (refIndex2 == -1? 0: 1);
        } else if (refIndex2 == -1) {
            return -1;
        }
        final int cmp = refIndex1 - refIndex2;
        if (cmp != 0) {
            return cmp;
        }
        return samRecord1.getMateAlignmentStart() - samRecord2.getMateAlignmentStart();
    }
}
