package au.edu.wehi.idsv;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;

import au.edu.wehi.idsv.bed.IntervalBed;

public class IntervalBedFilteringIterator<T extends DirectedEvidence> extends AbstractIterator<T> {
	private final IntervalBed filterOutRegions;
	private final Iterator<? extends T> it;
	private final int margin;
	public IntervalBedFilteringIterator(IntervalBed filterOutRegions, Iterator<? extends T> it, int margin) {
		if (filterOutRegions == null) throw new NullPointerException();
		if (it == null) throw new NullPointerException();
		if (margin < 0) throw new IllegalArgumentException();
		this.filterOutRegions = filterOutRegions;
		this.it = it;
		this.margin = margin;
	}
	@Override
	protected T computeNext() {
		while (it.hasNext()) {
			T r = it.next();
			BreakendSummary bs = r.getBreakendSummary();
			if (!filterOutRegions.overlaps(bs.referenceIndex, bs.start - margin, bs.end + margin)) {
				if (bs instanceof BreakpointSummary) {
					BreakpointSummary bp = (BreakpointSummary)bs;
					if (!filterOutRegions.overlaps(bp.referenceIndex2, bp.start2 - margin, bp.end2 + margin)) {
						return r;
					}
				} else {
					return r;
				}
			}
		}
		return endOfData();
	}
}
