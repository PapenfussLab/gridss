package scambly;

import au.edu.wehi.idsv.debruijn.PackedSequence;

import java.util.ArrayList;
import java.util.List;

public class ReadOverlap {
    Read left;
    Read right;
    /**
     * Number of bases the right read is offset from the start of the left read
     */
    int rightOffset;
    int overlapLength;
    int matchingBases;
    public ReadOverlap(Read left, Read right, int rightOffset) {
        this.left = left;
        this.right = right;
        this.rightOffset = rightOffset;
        this.overlapLength = PackedSequence.overlapLength(left.seq, right.seq, rightOffset);
        this.matchingBases = PackedSequence.overlapMatches(left.seq, right.seq, rightOffset);
    }
    public static List<ReadOverlap> findOverlaps(Read r1, Read r2) {
        List<ReadOverlap> list = new ArrayList<>();
        int r1len = r1.readLength();
        int r2len = r2.readLength();
        // [r1s, r1e]
        // [r2s, r2e]
        // r1s
        // r1               |---------|
        // r2  |-----|
        long max1OffsetIgnoringOverlap = r1.endFirstPosition - r2.startFirstPosition;
        // r1 |---------|
        // r2                |-----|
        long min1OffsetIgnoringOverlap = r1.startFirstPosition - r2.endFirstPosition;
        // TODO: restrict interval to bounds where the reads actually overlap
        long max2Before1 = r2.endFirstPosition - r1.startFirstPosition;

        long max1Offset = max1OffsetIgnoringOverlap;
        long min1Offset = min1OffsetIgnoringOverlap;
        for (int i = (int) min1Offset; i <= max1Offset; i++) {
            if (i < 0) {
                list.add(new ReadOverlap(r1, r2, -i));
            } else {
                list.add(new ReadOverlap(r2, r1, i));
            }
        }
        return list;
    }
}
