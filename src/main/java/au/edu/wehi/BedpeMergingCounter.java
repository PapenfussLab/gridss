package au.edu.wehi;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.util.IntervalUtil;
import htsjdk.samtools.util.Log;
import org.apache.commons.math3.util.Pair;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

public class BedpeMergingCounter {
    private static final Log log = Log.getInstance(BedpeMergingCounter.class);
    private static final BreakpointSummary SENTINEL = new BreakpointSummary(Integer.MAX_VALUE, BreakendDirection.Forward, 1, Integer.MAX_VALUE, BreakendDirection.Forward, 1);
    private static final int MERGE_MARGIN = 8;
    private final SortedMap<BreakpointSummary, Integer> active = new TreeMap<>(BreakpointSummary.ByStartStart2EndEnd2Direction);
    public List<Pair<BreakpointSummary, Integer>> process(BreakpointSummary bp) throws IOException {
        List<Pair<BreakpointSummary, Integer>> flushed = flushInactive(bp);
        if (!bp.isHighBreakend()) {
            process(bp, 1);
        }
        return flushed;
    }
    public List<Pair<BreakpointSummary, Integer>> finish() {
        return flushInactive(SENTINEL);
    }
    private void process(BreakpointSummary bp, int weight) throws IOException {
        // check if we already have this one
        Integer count = active.get(bp);
        if (count != null) {
            active.put(bp, count + weight);
            return;
        }
        for (BreakpointSummary existing : active.keySet()) {
            if (existing.direction == bp.direction && existing.direction2 == bp.direction2 &&
                    existing.referenceIndex == bp.referenceIndex && existing.referenceIndex2 == bp.referenceIndex2) {
                // we are contained in another interval
                if (existing.start <= bp.start && existing.start2 <= bp.start2 &&
                        existing.end >= bp.end && existing.end2 >= bp.end2) {
                    active.put(existing, active.get(existing) + weight);
                    return;
                }
                // we contain an existing interval
                if (existing.start >= bp.start && existing.start2 >= bp.start2 &&
                        existing.end <= bp.end && existing.end2 <= bp.end2) {
                    process(bp, active.remove(existing) + weight);
                    return;
                }
            }
        }
        // check if we can merge into an adjacent interval
        for (BreakpointSummary existing : active.keySet()) {
            if (existing.direction == bp.direction && existing.direction2 == bp.direction2 &&
                    existing.referenceIndex == bp.referenceIndex && existing.referenceIndex2 == bp.referenceIndex2) {
                if (existing.start == bp.start && existing.end == bp.end &&
                        IntervalUtil.overlapsClosed(existing.start2 - 1, existing.end2 + 1, bp.start2, bp.end2)) {
                    // merge intervals
                    BreakpointSummary merged = new BreakpointSummary(
                            existing.referenceIndex, existing.direction, existing.nominal, existing.start, existing.end,
                            existing.referenceIndex2, existing.direction2, existing.nominal2, Math.min(existing.start2, bp.start2), Math.max(existing.end2, bp.end2));
                    count = active.remove(existing);
                    process(merged, count + weight);
                    return;
                }
                if (existing.start2 == bp.start2 && existing.end2 == bp.end2 &&
                        IntervalUtil.overlapsClosed(existing.start - 1, existing.end + 1, bp.start, bp.end)) {
                    // merge intervals
                    BreakpointSummary merged = new BreakpointSummary(
                            existing.referenceIndex, existing.direction, existing.nominal, Math.min(existing.start, bp.start), Math.max(existing.end, bp.end),
                            existing.referenceIndex2, existing.direction2, existing.nominal2, existing.start2, existing.end2);
                    count = active.remove(existing);
                    process(merged, count + weight);
                    return;
                }
            }
        }
        // we're a new interval
        active.put(bp, weight);
    }
    private List<Pair<BreakpointSummary, Integer>> flushInactive(BreakpointSummary bp) {
        List<Pair<BreakpointSummary, Integer>> result = new ArrayList<>();
        while (!active.isEmpty()) {
            BreakpointSummary head = active.firstKey();
            if (head.referenceIndex < bp.referenceIndex || head.end < bp.start - MERGE_MARGIN) {
                int count = active.remove(head);
                result.add(Pair.create(head, count));
            } else {
                break;
            }
        }
        return result;
    }
}
