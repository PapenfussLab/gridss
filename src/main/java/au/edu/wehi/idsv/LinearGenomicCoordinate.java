package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

public interface LinearGenomicCoordinate {

	long getLinearCoordinate(int referenceIndex, long pos);

	long getLinearCoordinate(String chr, long pos);

	long getStartLinearCoordinate(BreakendSummary bp);

	long getStartLinearCoordinate(SAMRecord r);

	long getEndLinearCoordinate(BreakendSummary bp);

	long getEndLinearCoordinate(SAMRecord r);

	int getReferenceIndex(long linearCoordinate);

	int getReferencePosition(long linearCoordinate);

	String encodedToString(long linearCoordinate);

	String encodedIntervalToString(long startLinearCoordinate, long endLinearCoordinate);

	SAMSequenceDictionary getDictionary();
}