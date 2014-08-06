package au.edu.wehi.idsv.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

/**
 * Comparator for sorting SAMRecords by alignment coordinate end position
 */
public class SAMRecordEndCoordinateComparator extends SAMRecordCoordinateComparator {
	@Override
	public int fileOrderCompare(SAMRecord samRecord1, SAMRecord samRecord2) {
		final int refIndex1 = samRecord1.getReferenceIndex();
        final int refIndex2 = samRecord2.getReferenceIndex();
        if (refIndex1 == -1) {
            return (refIndex2 == -1? 0: 1);
        } else if (refIndex2 == -1) {
            return -1;
        }
        final int cmp = refIndex1 - refIndex2;
        if (cmp != 0) {
            return cmp;
        }
        return samRecord1.getAlignmentEnd() - samRecord2.getAlignmentEnd();
	}
}
