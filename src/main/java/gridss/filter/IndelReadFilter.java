package gridss.filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Filters out reads that do not contain any indels
 * @author cameron.d
 *
 */
public class IndelReadFilter implements SamRecordFilter {
	private final int minIndelLength;
	/**
	 * 
	 * @param minIndelLength minimum indel CIGAR element size
	 */
	public IndelReadFilter(int minIndelLength) {
		this.minIndelLength = minIndelLength;
	}
	public IndelReadFilter() {
		this(1);
	}
	@Override
	public boolean filterOut(SAMRecord record) {
		return record.getReadUnmappedFlag()
				|| record.getCigar() == null
				|| !record.getCigar().getCigarElements().stream().anyMatch(
						ce -> ce.getOperator().isIndelOrSkippedRegion() && ce.getLength() >= minIndelLength);
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		return filterOut(first) && filterOut(second);
	}

}
