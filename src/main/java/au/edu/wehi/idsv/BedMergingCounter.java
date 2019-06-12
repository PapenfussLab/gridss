package au.edu.wehi.idsv;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import org.apache.commons.math3.util.Pair;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

public class BedMergingCounter {
    private static final BreakendSummary SENTINEL = new BreakendSummary(Integer.MAX_VALUE, BreakendDirection.Forward, 1);
    public static final Ordering<BreakendSummary> ByStartEndDirection = new Ordering<BreakendSummary>() {
        public int compare(BreakendSummary o1, BreakendSummary o2) {
            return ComparisonChain.start()
                    .compare(o1.referenceIndex, o2.referenceIndex)
                    .compare(o1.start, o2.start)
                    .compare(o1.end, o2.end)
                    .compare(o1.direction, o2.direction)
                    .result();
        }};
    private final SortedMap<BreakendSummary, Integer> active = new TreeMap<>(ByStartEndDirection);
    private final boolean merge;

    public BedMergingCounter(boolean merge) {
        this.merge = merge;

    }
    public List<Pair<BreakendSummary, Integer>> process(BreakendSummary be) throws IOException {
        List<Pair<BreakendSummary, Integer>> flushed = flushInactive(be);
        process(be, 1);
        return flushed;
    }
    public List<Pair<BreakendSummary, Integer>> finish() {
        return flushInactive(SENTINEL);
    }
    private void process(BreakendSummary be, int weight) throws IOException {
        Integer existingCount = active.get(be);
        if (existingCount != null) {
            active.put(be, existingCount + weight);
            return;
        }
        if (merge) {
            BreakendSummary mergeTarget = new BreakendSummary(be.referenceIndex, be.direction, be.start - 1, be.start - 1, be.end + 1);
            for (BreakendSummary key : active.keySet()) {
                if (key.overlaps(mergeTarget)) {
                    BreakendSummary merged = new BreakendSummary(be.referenceIndex, be.direction,
                            key.start, Math.min(key.start, be.start), Math.max(key.end, be.end));
                    existingCount = active.remove(key);
                    process(merged, existingCount + weight);
                    return;
                }
            }
        }
        // new record
        active.put(be, weight);
    }
    private List<Pair<BreakendSummary, Integer>> flushInactive(BreakendSummary be) {
        List<Pair<BreakendSummary, Integer>> result = new ArrayList<>();
        while (!active.isEmpty()) {
            BreakendSummary head = active.firstKey();
            if (head.referenceIndex < be.referenceIndex || head.end < be.start - 1) {
                int count = active.remove(head);
                result.add(Pair.create(head, count));
            } else {
                break;
            }
        }
        return result;
    }
}
