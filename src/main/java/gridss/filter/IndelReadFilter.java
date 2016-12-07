package gridss.filter;

import java.util.List;

import htsjdk.samtools.CigarElement;
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
				|| !containsLargeIndel(record.getCigar().getCigarElements(), minIndelLength);
	}
	private static boolean containsLargeIndel(List<CigarElement> cigar, int minsize) {
		for (CigarElement ce : cigar) {
			if (ce.getLength() >= minsize && ce.getOperator().isIndelOrSkippedRegion()) {
				return true;
			}
		}
		return false;
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		return filterOut(first) && filterOut(second);
	}

}
