package au.edu.wehi.idsv;

import au.edu.wehi.idsv.debruijn.ReadErrorCorrector;
import com.google.common.collect.Iterables;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;

public class DirectedEvidenceErrorCorrectingIterator implements Iterator<DirectedEvidence> {
    private final LinearGenomicCoordinate linear;
    private final int k;
    private final float kmerErrorCorrectionMultiple;
    private final int bucketSize;
    private final Iterator<DirectedEvidence> in;
    private Deque<DirectedEvidence> lastBucket = new ArrayDeque<>();
    private Deque<DirectedEvidence> currentBucket = new ArrayDeque<>();
    private long currentBucketStart;

    public DirectedEvidenceErrorCorrectingIterator(
            LinearGenomicCoordinate linear,
            Iterator<DirectedEvidence> in,
            int minConcordantFragmentSize,
            int maxConcordantFragmentSize,
            int maxReadLength,
            int maxMappedReadLength,
            SAMEvidenceSource.EvidenceSortOrder eso,
            float kmerErrorCorrectionMultiple,
            int k) {
        if (eso != SAMEvidenceSource.EvidenceSortOrder.SAMRecordStartPosition) throw new IllegalArgumentException("NYI");
        int maxErrorCorrectSamRecordStartDelta = maxConcordantFragmentSize - minConcordantFragmentSize + 2 * maxReadLength - k;
        this.bucketSize = Math.max(maxErrorCorrectSamRecordStartDelta, 2 * maxMappedReadLength); 
        this.linear = linear;
        this.k = k;
        this.kmerErrorCorrectionMultiple = kmerErrorCorrectionMultiple;
        this.in = in;
        fillCurrentBucket();
    }

    @Override
    public boolean hasNext() {
        ensureLastBucket();
        return !lastBucket.isEmpty();
    }

    private void ensureLastBucket() {
        if (lastBucket.size() > 0) return;
        lastBucket = currentBucket;
        currentBucket = new ArrayDeque<>();
        fillCurrentBucket();
        // Error correct both buckets as since RPs could be spread across
        // multiple buckets due to differences in actual fragment size
        ReadErrorCorrector.errorCorrect(k, kmerErrorCorrectionMultiple, Iterables.concat(lastBucket, currentBucket));
    }

    private void fillCurrentBucket() {
        while (!addToCurrentBucket()) ;
    }

    /**
     * Adds to the current bucket
     * @return true if the current bucket is now full after adding, false otherwise
     */
    private boolean addToCurrentBucket() {
        if (!in.hasNext()) return true;
        DirectedEvidence de = in.next();
        SAMRecord underlying = de.getUnderlyingSAMRecord();
        if (currentBucket.isEmpty()) {
            currentBucketStart = linear.getStartLinearCoordinate(underlying);
        }
        currentBucket.add(de);
        return linear.getStartLinearCoordinate(underlying) - bucketSize > currentBucketStart;
    }

    @Override
    public DirectedEvidence next() {
        return lastBucket.removeFirst();
    }
}
