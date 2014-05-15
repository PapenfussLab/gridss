package au.edu.wehi.socrates;

import java.util.Comparator;

/**
 * Orders evidence by start position, then end position
 * If both records are of {@link DirectedBreakpoint}, the ordering sorts on the 
 * both start positions before comparing end position.
 *
 */
public class DirectedEvidenceCoordinateIntervalComparator implements Comparator<DirectedEvidence> {
	@Override
	public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
		BreakendSummary loc0 = arg0.getBreakendSummary();
		BreakendSummary loc1 = arg1.getBreakendSummary();
		int cmp = ComparatorUtil.compare(loc0.referenceIndex, loc1.referenceIndex);
		if (cmp == 0) cmp = ComparatorUtil.compare(loc0.start, loc1.start);
		if (cmp == 0 && loc0 instanceof BreakpointSummary && loc1 instanceof BreakpointSummary) {
			BreakpointSummary bi0 = (BreakpointSummary)loc0;
			BreakpointSummary bi1 = (BreakpointSummary)loc1;
			cmp = ComparatorUtil.compare(bi0.referenceIndex2, bi1.referenceIndex2);
			if (cmp == 0) cmp = ComparatorUtil.compare(bi0.start2, bi1.start2);
		}
		if (cmp == 0) cmp = ComparatorUtil.compare(loc0.end, loc1.end);
		if (cmp == 0 && loc0 instanceof BreakpointSummary && loc1 instanceof BreakpointSummary) {
			BreakpointSummary bi0 = (BreakpointSummary)loc0;
			BreakpointSummary bi1 = (BreakpointSummary)loc1;
			cmp = ComparatorUtil.compare(bi0.end2, bi1.end2);
		}
		return cmp;
	}
}
