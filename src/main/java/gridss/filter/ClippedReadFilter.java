package gridss.filter;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Filters out reads that are not hard or soft clipped
 * @author cameron.d
 *
 */
public class ClippedReadFilter implements SamRecordFilter {
	private final int minClipLength;
	/**
	 * 
	 * @param minClipLength minimum number of total clipped bases
	 */
	public ClippedReadFilter(int minClipLength) {
		this.minClipLength = minClipLength;
	}
	public ClippedReadFilter() {
		this(1);
	}
	@Override
	public boolean filterOut(SAMRecord record) {
		return record.getReadUnmappedFlag()
				|| record.getCigar() == null
				|| SAMRecordUtil.getEndClipLength(record) + SAMRecordUtil.getStartClipLength(record) < minClipLength;
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		return filterOut(first) && filterOut(second);
	}

}
