package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import java.util.Arrays;
import java.util.Objects;
import java.util.Set;

/**
 * Integer closed interval set
 */
public class IntegerIntervalSet {
    private static int DEFAULT_SIZE = 8;
    private int[] backing;
    private int count = 0;
    public IntegerIntervalSet() {
        backing = new int[DEFAULT_SIZE];
    }
    public IntegerIntervalSet(int[] pairedArray) {
        if ((pairedArray.length & 1) == 1) throw new IllegalArgumentException("Array does not contains integer pairs");
        for (int i = 1; i < pairedArray.length; i++) {
            if (pairedArray[i-1] > pairedArray[i]) {
                throw new IllegalArgumentException("Interval pairs are not sorted in ascending order");
            }
        }
        this.backing = pairedArray;
        this.count = pairedArray.length;
    }

    public IntegerIntervalSet(int size) {
        backing = new int[2 * size];
    }

    public void add(int low, int high) {
        int startOffset = 0;
        while (startOffset < count && intervalEnd(startOffset) < low) {
            startOffset++;
        }
        if (startOffset == count) {
            // after last interval
            append(low, high);
        } else if (high < intervalStart(startOffset)) {
            // between intervals
            ensureSize(count + 1);
            System.arraycopy(backing, 2 * startOffset, backing, 2 * startOffset + 2, (count - startOffset) * 2);
            backing[2 * startOffset] = low;
            backing[2 * startOffset + 1] = high;
            count++;
        } else {
            // overlap an interval
            backing[2 * startOffset] = Math.min(low, intervalStart(startOffset));
            int endOffset = startOffset;
            while (endOffset < count && intervalStart(endOffset) <= high) {
                endOffset++;
            }
            // merge overlapping intervals
            endOffset--;
            assert(startOffset <= endOffset);
            backing[2 * startOffset + 1] = Math.max(high, intervalEnd(endOffset));
            int shiftBy = endOffset - startOffset;
            int recordsToShift = count - endOffset - 1;
            if (shiftBy > 0 && recordsToShift > 0) {
                // shuffle all the records after our final merged interval into their new location
                System.arraycopy(backing, 2 * (endOffset + 1), backing, 2 * (startOffset + 1), 2 * recordsToShift);
            }
            count -= shiftBy;
        }
    }

    private void ensureSize(int size) {
        if ((backing.length >> 1) < size) {
            int[] arr = new int[2 * backing.length];
            System.arraycopy(backing, 0, arr, 0, 2 * count);
            backing = arr;
        }
    }

    public void append(int low, int high) {
        if (count > 0 && backing[2 * (count - 1) + 1] > low) {
            throw new IllegalArgumentException(String.format("Attempting to append [%d, %d] to %s", low, high, toString()));
        }
        ensureSize(count + 1);
        backing[2 * count] = low;
        backing[2 *count + 1] = high;
        count++;
    }

    public boolean encloses(int low, int high) {
        for (int i = 0; i < count; i++) {
            if (intervalStart(i) <= low && intervalEnd(i) >= high) {
                return true;
            }
        }
        return false;
    }
    public int intervalStart(int i) {
        return backing[2 * i];
    }
    public int intervalEnd(int i) {
        return backing[2 * i + 1];
    }
    public int size() {
        return count;
    }
    public boolean isEmpty() {
        return count == 0;
    }
    public static IntegerIntervalSet intersect(IntegerIntervalSet set1, IntegerIntervalSet set2) {
        IntegerIntervalSet out = new IntegerIntervalSet();
        int offset1 = 0;
        int offset2 = 0;
        while(offset1 < set1.count && offset2 < set2.count) {
            int start1 = set1.intervalStart(offset1);
            int end1 = set1.intervalEnd(offset1);
            int start2 = set2.intervalStart(offset2);
            int end2 = set2.intervalEnd(offset2);
            if (end1 < start2) {
                offset1++;
            } else if (end2 < start1) {
                offset2++;
            } else {
                int start = Math.max(start1, start2);
                int end = Math.min(end1, end2);
                out.append(start, end);
                if (end1 < end2) {
                    offset1++;
                } else {
                    offset2++;
                }
            }
        }
        return out;
    }

    public Set<Range<Integer>> asRanges() {
        RangeSet<Integer> rs = TreeRangeSet.create();
        for (int i = 0; i < size(); i++) {
            rs.add(Range.closed(intervalStart(i), intervalEnd(i)));
        }
        return rs.asRanges();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("{");
        for (int i = 0; i < count; i++) {
            sb.append('[');
            sb.append(intervalStart(i));
            sb.append(',');
            sb.append(intervalEnd(i));
            sb.append(']');
            sb.append(", ");
        }
        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof IntegerIntervalSet)) return false;
        IntegerIntervalSet that = (IntegerIntervalSet) o;
        return count == that.count &&
                arraysEqualUntil(backing, that.backing, count);
    }
    private static boolean arraysEqualUntil(int[] a, int[] b, int n) {
        for (int i = 0; i < n; i++) {
            if (a[n] != b[n]) {
                return false;
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        int result = count;
        result = 31 * result + Arrays.hashCode(backing);
        return result;
    }
}
