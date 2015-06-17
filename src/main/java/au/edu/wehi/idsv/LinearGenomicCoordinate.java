package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public interface LinearGenomicCoordinate {

	public abstract long getLinearCoordinate(int referenceIndex, long pos);

	public abstract long getLinearCoordinate(String chr, long pos);

	public abstract long getStartLinearCoordinate(BreakendSummary bp);

	public abstract long getStartLinearCoordinate(SAMRecord r);

	public abstract long getEndLinearCoordinate(BreakendSummary bp);

	public abstract long getEndLinearCoordinate(SAMRecord r);

	public abstract int getReferenceIndex(long linearCoordinate);

	public abstract int getReferencePosition(long linearCoordinate);

	public abstract String encodedToString(long linearCoordinate);

	public abstract String encodedIntervalToString(long startLinearCoordinate,
			long endLinearCoordinate);

}