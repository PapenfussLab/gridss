package au.edu.wehi.socrates;

import java.util.Comparator;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordComparator;
import net.sf.samtools.SAMRecordCoordinateComparator;

/**
 * Orders evidence by start position, then end position
 * If both records are of {@link DirectedBreakpoint}, the ordering sorts on the 
 * both start positions before comparing end position.
 *
 */
public class DirectedEvidenceCoordinateIntervalComparator implements Comparator<DirectedEvidence> {
	@Override
	public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
		BreakpointLocation loc0 = arg0.getBreakpointLocation();
		BreakpointLocation loc1 = arg1.getBreakpointLocation();
		int cmp = ComparatorUtil.integerCompare(loc0.referenceIndex, loc1.referenceIndex);
		if (cmp == 0) cmp = ComparatorUtil.integerCompare(loc0.start, loc1.start);
		if (cmp == 0 && loc0 instanceof BreakpointInterval && loc1 instanceof BreakpointInterval) {
			BreakpointInterval bi0 = (BreakpointInterval)loc0;
			BreakpointInterval bi1 = (BreakpointInterval)loc1;
			cmp = ComparatorUtil.integerCompare(bi0.referenceIndex2, bi1.referenceIndex2);
			if (cmp == 0) cmp = ComparatorUtil.integerCompare(bi0.start2, bi1.start2);
		}
		if (cmp == 0) cmp = ComparatorUtil.integerCompare(loc0.end, loc1.end);
		if (cmp == 0 && loc0 instanceof BreakpointInterval && loc1 instanceof BreakpointInterval) {
			BreakpointInterval bi0 = (BreakpointInterval)loc0;
			BreakpointInterval bi1 = (BreakpointInterval)loc1;
			cmp = ComparatorUtil.integerCompare(bi0.end2, bi1.end2);
		}
		return cmp;
	}
}
