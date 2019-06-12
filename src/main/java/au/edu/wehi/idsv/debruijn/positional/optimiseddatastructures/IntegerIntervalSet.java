package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

/**
 * Integer closed interval set
 */
public class IntegerIntervalSet {
    private static boolean SANITY_CHECK_ARRAY_BOUNDS = true;
    private static int DEFAULT_SIZE = SANITY_CHECK_ARRAY_BOUNDS ? 0 : 8;
    private int[] backing = new int[DEFAULT_SIZE];
    private int count = 0;
    public IntegerIntervalSet() {
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
        if (count < size) {
            resizeArray(SANITY_CHECK_ARRAY_BOUNDS ? (2 * size) : (2 * count));
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

    private void resizeArray(int targetSize) {
        int[] arr = new int[targetSize];
        System.arraycopy(backing, 0, arr, 0, 2 * count);
        backing = arr;
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
    public IntegerIntervalSet intersect(IntegerIntervalSet set) {
        throw new RuntimeException("NYI");
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
}
