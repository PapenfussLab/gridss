package au.edu.wehi.idsv.util;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

public class RangeUtil {
	public static RangeSet<Integer> intersect(RangeSet<Integer> rs1, RangeSet<Integer> rs2) {
		RangeSet<Integer> rs = TreeRangeSet.create();
		for (Range<Integer> r : rs1.asRanges()) {
			rs.addAll(rs2.subRangeSet(r));
		}
		return rs;
	}
}
