package au.edu.wehi.idsv.util;

import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Map.Entry;

import static org.junit.Assert.assertEquals;

public class RangeUtilTest {
	@Test
	public void addWhereBest_should_insert_when_best() {
		RangeMap<Integer, Integer> rm = TreeRangeMap.create();
		RangeUtil.addWhereBest(rm, Range.closedOpen(0, 10), 1, Comparator.naturalOrder());
		RangeUtil.addWhereBest(rm, Range.closedOpen(1, 5), 1, Comparator.naturalOrder());
		RangeUtil.addWhereBest(rm, Range.closedOpen(4, 7), 4, Comparator.naturalOrder());
		RangeUtil.addWhereBest(rm, Range.closedOpen(6, 10), 3, Comparator.naturalOrder());
		RangeUtil.addWhereBest(rm, Range.closedOpen(8, 9), 5, Comparator.naturalOrder());
		ArrayList<Entry<Range<Integer>, Integer>> result = Lists.newArrayList(rm.asMapOfRanges().entrySet());
		assertEquals(Range.closedOpen(0, 4), result.get(0).getKey());
		assertEquals(1, (int)result.get(0).getValue());
		assertEquals(Range.closedOpen(4, 7), result.get(1).getKey());
		assertEquals(4, (int)result.get(1).getValue());
		assertEquals(Range.closedOpen(7, 8), result.get(2).getKey());
		assertEquals(3, (int)result.get(2).getValue());
		assertEquals(Range.closedOpen(8, 9), result.get(3).getKey());
		assertEquals(5, (int)result.get(3).getValue());
		assertEquals(Range.closedOpen(9, 10), result.get(4).getKey());
		assertEquals(3, (int)result.get(4).getValue());
	}
	// Just a sanity check test case to ensure that my reading of the guava docs is correct
	@Test
	public void remove_should_split_ranges() {
		RangeMap<Integer, Integer> rm = TreeRangeMap.create();
		rm.put(Range.closedOpen(0, 10), 1);
		rm.remove(Range.lessThan(5));
		assertEquals(Range.closedOpen(5, 10), rm.getEntry(5).getKey());
	}
}
