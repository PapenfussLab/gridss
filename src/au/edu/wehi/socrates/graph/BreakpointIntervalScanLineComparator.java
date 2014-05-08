package au.edu.wehi.socrates.graph;

import java.util.Comparator;

import au.edu.wehi.socrates.BreakpointInterval;
import au.edu.wehi.socrates.BreakpointLocation;
import au.edu.wehi.socrates.ComparatorUtil;

/**
 * Orders breakpoints according to a scan-line traversal of 2D
 * projection of the linearised coordinates of the breakpoints
 * @author Daniel Cameron
 *
 */
public class BreakpointIntervalScanLineComparator implements Comparator<BreakpointInterval> {
	@Override
	public int compare(BreakpointInterval loc0, BreakpointInterval loc1) {
		int cmp = ComparatorUtil.compare(loc0.referenceIndex, loc1.referenceIndex);
		if (cmp == 0) cmp = ComparatorUtil.compare(loc0.start, loc1.start);
		if (cmp == 0) cmp = ComparatorUtil.compare(loc0.referenceIndex2, loc1.referenceIndex2);
		if (cmp == 0) cmp = ComparatorUtil.compare(loc0.start2, loc1.start2);
		if (cmp == 0) cmp = ComparatorUtil.compare(loc0.end, loc1.end);
		if (cmp == 0) cmp = ComparatorUtil.compare(loc0.end2, loc1.end2);
		return cmp;
	}

}
