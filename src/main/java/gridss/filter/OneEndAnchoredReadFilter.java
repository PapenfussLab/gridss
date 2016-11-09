package gridss.filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Filters out reads that are not part of a read pair in which
 * only one read in the read pair is mapped.
 * @author Daniel Cameron
 *
 */
public class OneEndAnchoredReadFilter implements SamRecordFilter {

	@Override
	public boolean filterOut(SAMRecord record) {
		boolean isOEA = record.getReadPairedFlag() && (record.getReadUnmappedFlag() ^ record.getMateUnmappedFlag());
		return !isOEA;
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		// both have to be OEA
		return filterOut(first) || filterOut(second);
	}
}
