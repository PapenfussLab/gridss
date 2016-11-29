package gridss.filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Filters all read reads
 * @author Daniel Cameron
 *
 */
public class FixedFilter implements SamRecordFilter {
	private final boolean filter;
	public FixedFilter(boolean filterOut) {
		this.filter = filterOut;
	}
	@Override
	public boolean filterOut(SAMRecord record) {
		return filter;
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		return filter;
	}
}
