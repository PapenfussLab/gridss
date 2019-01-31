package au.edu.wehi.idsv;

import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class BedMergingCounterTest extends TestHelper {
    @Test
    public void shouldMerge() throws IOException {
        List<Pair<BreakendSummary, Integer>> result = new ArrayList<>();
        BedMergingCounter c = new BedMergingCounter(true);
        result.addAll(c.process(new BreakendSummary(0, FWD, 1)));
        result.addAll(c.process(new BreakendSummary(0, FWD, 2)));
        result.addAll(c.process(new BreakendSummary(0, FWD, 3)));
        result.addAll(c.process(new BreakendSummary(0, FWD, 5)));
        result.addAll(c.process(new BreakendSummary(1, FWD, 5)));
        result.addAll(c.finish());
        assertEquals(3, result.size());
        assertEquals(3, (int)result.get(0).getSecond());
        assertEquals(new BreakendSummary(0, FWD, 1, 1, 3), result.get(0).getFirst());
        assertEquals(1, (int)result.get(1).getSecond());
        assertEquals(1, (int)result.get(2).getSecond());
    }

}