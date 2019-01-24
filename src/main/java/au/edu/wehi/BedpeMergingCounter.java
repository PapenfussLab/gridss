package au.edu.wehi;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakpointSummary;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import htsjdk.samtools.util.Log;
import org.apache.commons.math3.util.Pair;

import java.io.IOException;
import java.util.*;

public class BedpeMergingCounter {
    private static final Log log = Log.getInstance(BedpeMergingCounter.class);
    private static final BreakpointSummary SENTINEL = new BreakpointSummary(Integer.MAX_VALUE, BreakendDirection.Forward, 1, Integer.MAX_VALUE, BreakendDirection.Forward, 1);
    private static final int MERGE_MARGIN = 2;
    private final SortedSet<BreakpointSummary> activeByEnd1 = new TreeSet<>(ByEndStartEnd2Start2Direction12);
    // lookup on the remote location since we're doing a linear traversal of the local thus it is uninformative
    private final SortedMap<BreakpointSummary, Integer> activeByEnd2 = new TreeMap<>(ByEnd2Start2EndIStartDirection12);
    // how far around the breakpoint position we need to check to ensure we find any potential overlaps
    private int maxWidth2 = 0;
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
        maxWidth2 = Math.max(maxWidth2, bp.end2 - bp.start2 + 1);
        if (activeByEnd1.contains(bp)) {
            activeByEnd2.put(bp, activeByEnd2.get(bp) + weight);
            return;
        }
        // check for overlaps or adjacencies
        BreakpointSummary lookup = new BreakpointSummary(
                bp.referenceIndex, bp.direction, bp.start - 1, bp.start - 1, bp.end + 1,
                bp.referenceIndex2, bp.direction2, bp.start2 - 1, bp.start2 - 1, bp.end2 + 1);

        // lookup the start of the potential overlap or adjacency interval
        // as we could potentially overlap anything that ends after our starting position
        SortedMap<BreakpointSummary, Integer> overlap2 = activeByEnd2.tailMap(new BreakpointSummary(
                bp.referenceIndex, bp.direction, bp.start - 1,
                bp.referenceIndex2, bp.direction2, bp.start2 - 1 - maxWidth2))
                .headMap(new BreakpointSummary(
                        bp.referenceIndex, bp.direction, bp.end + 1,
                        bp.referenceIndex2, bp.direction2, bp.end2 + 1 + maxWidth2));
        for (BreakpointSummary key : overlap2.keySet()) {
            if (key.overlaps(lookup)) {
                BreakpointSummary merged = new BreakpointSummary(
                        key.referenceIndex, key.direction, key.nominal, Math.min(key.start, bp.start), Math.max(key.end, bp.end),
                        key.referenceIndex2, key.direction2, key.nominal2, Math.min(key.start2, bp.start2), Math.max(key.end2, bp.end2));
                int existingWeight = activeByEnd2.remove(key);
                activeByEnd1.remove(key);
                process(merged, existingWeight + weight);
                return;
            }
        }
        // we're a new interval
        activeByEnd1.add(bp);
        activeByEnd2.put(bp, weight);
    }
    private List<Pair<BreakpointSummary, Integer>> flushInactive(BreakpointSummary bp) {
        List<Pair<BreakpointSummary, Integer>> result = new ArrayList<>();
        while (!activeByEnd1.isEmpty()) {
            BreakpointSummary head = activeByEnd1.first();
            if (head.referenceIndex < bp.referenceIndex || head.end < bp.start - MERGE_MARGIN) {
                activeByEnd1.remove(head);
                int count = activeByEnd2.remove(head);
                result.add(Pair.create(head, count));
            } else {
                break;
            }
        }
        return result;
    }
    // ordering both ignore nominal since we don't care about it
    public static Ordering<BreakpointSummary> ByEndStartEnd2Start2Direction12 = new Ordering<BreakpointSummary>() {
        public int compare(BreakpointSummary o1, BreakpointSummary o2) {
            return ComparisonChain.start()
                    .compare(o1.referenceIndex, o2.referenceIndex)
                    .compare(o1.end, o2.end)
                    .compare(o1.start, o2.start)
                    .compare(o1.referenceIndex2, o2.referenceIndex2)
                    .compare(o1.end2, o2.end2)
                    .compare(o1.start2, o2.start2)
                    .compare(o1.direction, o2.direction)
                    .compare(o1.direction2, o2.direction2)
                    .result();
        }
    };
    public static Ordering<BreakpointSummary> ByEnd2Start2EndIStartDirection12 = new Ordering<BreakpointSummary>() {
        public int compare(BreakpointSummary o1, BreakpointSummary o2) {
            return ComparisonChain.start()
                    .compare(o1.referenceIndex2, o2.referenceIndex2)
                    .compare(o1.end2, o2.end2)
                    .compare(o1.start2, o2.start2)
                    .compare(o1.referenceIndex, o2.referenceIndex)
                    .compare(o1.end, o2.end)
                    .compare(o1.start, o2.start)
                    .compare(o1.direction, o2.direction)
                    .compare(o1.direction2, o2.direction2)
                    .result();
        }
    };
}
