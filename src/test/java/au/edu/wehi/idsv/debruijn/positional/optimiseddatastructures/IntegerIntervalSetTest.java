package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.util.RangeUtil;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class IntegerIntervalSetTest extends TestHelper {
    private final List<Range<Integer>> ALL_RANGES_10 = allRanges(0, 10);
    private final List<RangeSet<Integer>> ALL_PAIRS_10 = allRangePairs(0, 10);
    public List<Range<Integer>> allRanges(int low, int high) {
        List<Range<Integer>> list = new ArrayList<>((high + low) / 2);
        for (int i = low; i <= high; i++) {
            for (int j = i; j <= high; j++) {
                list.add(Range.closed(i, j));
            }
        }
        return list;
    }
    public List<RangeSet<Integer>> allRangePairs(int low, int high) {
        List<RangeSet<Integer>> list = new ArrayList<>();
        for (int i = low; i <= high; i++) {
            for (int j = i + 1; j <= high; j++) {
                for (int k = j + 1; k <= high; k++) {
                    for (int l = k + 1; l <= high; l++) {
                        list.add(RS(i, j, k, l));
                    }
                }
            }
        }
        return list;
    }
    @Test
    public void should_mimic_rangeset() {
        for (Range<Integer> r1 : ALL_RANGES_10) {
            for (Range<Integer> r2 : ALL_RANGES_10) {
                assertMatching(ImmutableList.of(r1, r2), true);
                for (Range<Integer> r3 : ALL_RANGES_10) {
                    for (Range<Integer> r4 : ALL_RANGES_10) {
                        assertMatching(ImmutableList.of(r1, r2, r3, r4), false);
                    }
                }
            }
        }
    }
    @Test
    public void intersect_should_mimic_rangeset() {
        for (RangeSet<Integer> r1 : ALL_PAIRS_10) {
            for (RangeSet<Integer> r2 : ALL_PAIRS_10) {
                testIntersect(r1, r2);
            }
            for (Range<Integer> r2 : ALL_RANGES_10) {
                testIntersect(r1, TreeRangeSet.create(ImmutableList.of(r2)));
                testIntersect(TreeRangeSet.create(ImmutableList.of(r2)), r1);
            }
        }
    }

    /**
     * Matches RangeSet after adding the given intervals in the given order
     * @param list
     */
    private void assertMatching(List<Range<Integer>> list, boolean assertEncloses) {
        RangeSet<Integer> rs = TreeRangeSet.create();
        IntegerIntervalSet iis = new IntegerIntervalSet();
        for (Range<Integer> r : list) {
            rs.add(r);
            iis.add(r.lowerEndpoint(), r.upperEndpoint());
            assertMatching(rs, iis, assertEncloses);
        }
    }

    private void assertMatching(RangeSet<Integer> rs, IntegerIntervalSet iis, boolean assertEncloses) {
        assertEquals(rs.asRanges().size(), iis.size());
        int i = 0;
        for (Range<Integer> r : rs.asRanges()) {
            assertEquals((int)r.lowerEndpoint(), iis.intervalStart(i));
            assertEquals((int)r.upperEndpoint(), iis.intervalEnd(i));
            i++;
        }
        if (assertEncloses) {
            for (Range<Integer> r : ALL_RANGES_10) {
                assertEquals(rs.encloses(r), iis.encloses(r.lowerEndpoint(), r.upperEndpoint()));
            }
        }
    }

    private void testIntersect(RangeSet<Integer> set1, RangeSet<Integer> set2) {
        RangeSet<Integer> rs = RangeUtil.intersect(set1, set2);
        IntegerIntervalSet iis = IntegerIntervalSet.intersect(IIS(set1), IIS(set2));
        assertMatching(rs, iis, false);
    }
    private Range<Integer> r(int low, int high) { return Range.closed(low, high);}
    @Test
    public void testIntervalMerge() {
        assertMatching(ImmutableList.of(
                r(0, 0),
                r(0,1),
                r(0,0)), true);
        assertMatching(ImmutableList.of(
                r(0, 0),
                r(1,1),
                r(0,1)), true);
        assertMatching(ImmutableList.of(
                r(0, 0),
                r(1,1),
                r(2,2),
                r(3,3),
                r(4,4),
                r(0, 4)), true);
        assertMatching(ImmutableList.of(
                r(0, 0),
                r(10,10),
                r(20,20),
                r(30,30),
                r(40,40),
                r(5, 35)), true);
        assertMatching(ImmutableList.of(
                r(0, 0),
                r(10,10),
                r(20,20),
                r(30,30),
                r(40,40),
                r(5, 40)), true);
    }
    @Test
    public void intersect() {
        testIntersect(RS(0, 1, 2, 3, 4, 5), RS(-1, 6));
        testIntersect(RS(0, 1, 2, 3, 4, 5), RS(0, 6));
        testIntersect(RS(0, 1, 2, 3, 4, 5), RS(1, 6));
        testIntersect(RS(0, 1, 2, 3, 4, 5), RS(0, 5));
        testIntersect(RS(0, 1, 2, 3, 4, 5), RS(0, 4));
    }
}