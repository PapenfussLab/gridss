package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

/**
 * Compares SAMRecord by coordinate only. For consistency with SAM/BAM sort ordering,
 * records with no alignment reported are considered to be after all aligned records.
 */
public class SAMRecordCoordinateOnlyComparator extends SAMRecordCoordinateComparator {
    @Override
    public int compare(SAMRecord o1, SAMRecord o2) {
        return fileOrderCompare(o1, o2);
    }
}
