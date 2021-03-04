package gridss.filter;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Filters out reads that are not hard or soft clipped
 * @author cameron.d
 *
 */
public class ClippedReadFilter implements SamRecordFilter {
	private final int minClipLength;
	private final boolean includeSplitReads;
	/**
	 * 
	 * @param minClipLength minimum number of total clipped bases
	 */
	public ClippedReadFilter(int minClipLength, boolean includeSplitReads) {
		this.minClipLength = minClipLength;
		this.includeSplitReads = includeSplitReads;
	}
	public ClippedReadFilter() {
		this(1, false);
	}
	@Override
	public boolean filterOut(SAMRecord record) {
		return record.getReadUnmappedFlag()
				|| record.getCigar() == null
				|| (!includeSplitReads && record.hasAttribute(SAMTag.SA.name()))
				|| SAMRecordUtil.getEndClipLength(record) + SAMRecordUtil.getStartClipLength(record) < minClipLength;
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		return filterOut(first) && filterOut(second);
	}

}
