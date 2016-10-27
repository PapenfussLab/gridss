package gridss.filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Filters out reads that have a single alignment
 * @author Daniel Cameron
 *
 */
public class SplitReadFilter implements SamRecordFilter {

	@Override
	public boolean filterOut(SAMRecord record) {
		return record.getAttribute(SAMTag.SA.name()) == null;
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		return filterOut(first) && filterOut(second);
	}
}
