package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;


/**
 * Offsets the given linear genomic coordinate system by a fixed amount,
 * reducing the range valid values to an integer.
 * 
 * Throws RuntimeException whenever the range would fall outside the valid signed 32 bit integer range. 
 * 
 * @author Daniel Cameron
 *
 */
public class IntegerLinearGenomicCoordinate implements LinearGenomicCoordinate {
	private final LinearGenomicCoordinate underlying;
	private final long offset;
	public IntegerLinearGenomicCoordinate(LinearGenomicCoordinate underlying, long offset) {
		this.underlying = underlying;
		this.offset = offset;
	}
	private long removeOffset(long withOffset) {
		long withoutOffset = withOffset - offset; 
		if (withOffset > Integer.MAX_VALUE || withOffset < Integer.MIN_VALUE) {
			throw new RuntimeException(String.format("%s out of integer range given offset %s", underlying.encodedToString(withoutOffset), underlying.encodedToString(offset)));
		}
		return withoutOffset;
	}
	private long addOffset(long withoutOffset) {
		long withOffset = withoutOffset + offset;
		if (withOffset > Integer.MAX_VALUE || withOffset < Integer.MIN_VALUE) {
			throw new RuntimeException(String.format("%s out of integer range given offset %s", underlying.encodedToString(withoutOffset), underlying.encodedToString(offset)));
		}
		return withOffset;
	}
	@Override
	public long getLinearCoordinate(int referenceIndex, long pos) {
		return addOffset(underlying.getLinearCoordinate(referenceIndex, pos));
	}
	@Override
	public long getLinearCoordinate(String chr, long pos) {
		return addOffset(underlying.getLinearCoordinate(chr, pos));
	}
	@Override
	public long getStartLinearCoordinate(BreakendSummary bp) {
		return addOffset(underlying.getStartLinearCoordinate(bp));
	}
	@Override
	public long getStartLinearCoordinate(SAMRecord r) {
		return addOffset(underlying.getStartLinearCoordinate(r));
	}
	@Override
	public long getEndLinearCoordinate(BreakendSummary bp) {
		return addOffset(underlying.getEndLinearCoordinate(bp));
	}
	@Override
	public long getEndLinearCoordinate(SAMRecord r) {
		return addOffset(underlying.getEndLinearCoordinate(r));
	}
	@Override
	public int getReferenceIndex(long linearCoordinate) {
		return underlying.getReferenceIndex(removeOffset(linearCoordinate));
	}
	@Override
	public int getReferencePosition(long linearCoordinate) {
		return underlying.getReferencePosition(removeOffset(linearCoordinate));
	}
	@Override
	public String encodedToString(long linearCoordinate) {
		return underlying.encodedToString(removeOffset(linearCoordinate));
	}
	@Override
	public String encodedIntervalToString(long startLinearCoordinate, long endLinearCoordinate) {
		return underlying.encodedIntervalToString(removeOffset(startLinearCoordinate), removeOffset(endLinearCoordinate));
	}
}
