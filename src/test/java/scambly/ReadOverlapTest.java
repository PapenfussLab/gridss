package scambly;

import au.edu.wehi.idsv.TestHelper;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;

@Ignore("Unfinished prototype")
public class ReadOverlapTest extends TestHelper {
    public Read R(int length, int start, int end) {
        return new Read(getContext().getLinear(), Read(0, 1, String.format("%dM", length)), start, end);
    }
    @Test
    public void findOverlaps_should_report_all_overlaps() {
        List<ReadOverlap> overlaps = ReadOverlap.findOverlaps(
                R(50, 100, 200),
                R(50, 100, 200));
        assertEquals(2 * 49, overlaps.size());
    }
    @Test
    public void findOverlaps_should_restrict_overlaps_to_valid_positions() {
        assertEquals(0,ReadOverlap.findOverlaps(R(1, 1, 1), R(1, 2, 2)).size());
        assertEquals(1,ReadOverlap.findOverlaps(R(1, 1, 1), R(1, 1, 1)).size());
        assertEquals(1, ReadOverlap.findOverlaps(R(2, 0, 1), R(1, 1, 1)).size());
        assertEquals(2, ReadOverlap.findOverlaps(R(2, 0, 1), R(2, 1, 1)).size());
        assertEquals(3, ReadOverlap.findOverlaps(R(2, 0, 2), R(2, 1, 1)).size());
        assertEquals(3, ReadOverlap.findOverlaps(R(2, 0, 3), R(2, 1, 1)).size());
        assertEquals(1, ReadOverlap.findOverlaps(R(1, 0, 100), R(1, 1, 100)).size());
        assertEquals(2, ReadOverlap.findOverlaps(R(1, 0, 100), R(2, 1, 100)).size());
        assertEquals(3, ReadOverlap.findOverlaps(R(2, 0, 100), R(2, 1, 100)).size());
        assertEquals(4, ReadOverlap.findOverlaps(R(2, 0, 100), R(3, 1, 100)).size());
        assertEquals(5, ReadOverlap.findOverlaps(R(2, 0, 100), R(4, 1, 100)).size());
        assertEquals(6, ReadOverlap.findOverlaps(R(3, 0, 100), R(4, 1, 100)).size());
        assertEquals(7, ReadOverlap.findOverlaps(R(4, 0, 100), R(4, 1, 100)).size());
    }
}