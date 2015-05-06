package au.edu.wehi.idsv;

import java.util.Iterator;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.bed.IntervalBed;

import com.google.common.collect.AbstractIterator;

public class IntervalBedFilteringIterator extends AbstractIterator<SAMRecord> {
	private final IntervalBed filterOutRegions;
	private final Iterator<SAMRecord> it;
	private final int margin;
	public IntervalBedFilteringIterator(IntervalBed filterOutRegions, Iterator<SAMRecord> it, int margin) {
		this.filterOutRegions = filterOutRegions;
		this.it = it;
		this.margin = margin;
	}
	@Override
	protected SAMRecord computeNext() {
		while (it.hasNext()) {
			SAMRecord r = it.next();
			if (r.getReadUnmappedFlag() || !filterOutRegions.overlaps(r.getReferenceIndex(), r.getAlignmentStart() - margin, r.getAlignmentEnd() + margin)) {
				return r;
			}
		}
		return endOfData();
	}
}
