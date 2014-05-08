package au.edu.wehi.socrates.graph;

import java.util.Comparator;

import au.edu.wehi.socrates.ComparatorUtil;

/**
 * Orders according to endX, endY, startX, startY 
 * @author cameron.d
 *
 */
public class TrapezoidGraphNodeEndXYComparator implements Comparator<TrapezoidGraphNode> {
	@Override
	public int compare(TrapezoidGraphNode o1, TrapezoidGraphNode o2) {
		int cmp = ComparatorUtil.compare(o1.endX, o2.endX);
		if (cmp == 0) cmp = ComparatorUtil.compare(o1.endY, o2.endY);
		if (cmp == 0) cmp = ComparatorUtil.compare(o1.startX, o2.startX);
		if (cmp == 0) cmp = ComparatorUtil.compare(o1.startY, o2.startY);
		return cmp;
	}
}