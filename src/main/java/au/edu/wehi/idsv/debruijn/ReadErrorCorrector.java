package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.longs.Long2IntMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2LongMap;
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap;

import java.util.HashSet;
import java.util.Set;

public class ReadErrorCorrector {
    private static long SENTINEL_VALUE = Long.MIN_VALUE;
    private static int MAX_BASE_CORRECTIONS = 2;
    private static final Log log = Log.getInstance(ReadErrorCorrector.class);
    private final Long2IntMap kmerCounts = new Long2IntOpenHashMap();
    private final int k;
    private final float collapseMultiple;
    private Long2LongMap collapseLookup;

    public ReadErrorCorrector(int k, float collapseMultiple) {
        if (k > 31) throw new IllegalArgumentException("k cannot exceed 31");
        this.k = k;
        this.collapseMultiple = collapseMultiple;
    }

    public static void errorCorrect(int k, float collapseMultiple, Iterable<? extends DirectedEvidence> evidence) {
        ReadErrorCorrector ec = new ReadErrorCorrector(k, collapseMultiple);
        // need to deduplicate the underlying reads so we don't double count
        // kmers from reads with multiple evidence (e.g. multiple indels or SC on both ends)
        Set<SAMRecord> reads = new HashSet<>();
        for (DirectedEvidence de : evidence) {
            reads.add(de.getUnderlyingSAMRecord());
            if (de instanceof NonReferenceReadPair) {
                reads.add(((NonReferenceReadPair) de).getNonReferenceRead());
            }
        }
        reads.stream().forEach(r -> ec.countKmers(r));
        reads.stream().forEach(r -> ec.errorCorrect(r));
    }

    public void countKmers(SAMRecord r) {
        PackedSequence ps = new PackedSequence(r.getReadBases(), false, false);
        for (int i = 0; i < ps.length() - k + 1; i++) {
            long kmer = ps.getKmer(i, k);
            int count = kmerCounts.get(kmer);
            count++;
            kmerCounts.put(kmer, count);
        }
        collapseLookup = null;
    }
    public int errorCorrect(SAMRecord r) {
        if (r.getReadLength() < k) return 0;
        ensureCollapseLookup();
        PackedSequence ps = new PackedSequence(r.getReadBases(), false, false);
        int changes = error_correct(r, ps);
        if (changes > 0) {
            r.setReadBases(ps.getBytes(0, r.getReadBases().length));
        }
        return changes;
    }
    public int error_correct(SAMRecord r, PackedSequence ps) {
        int changes = 0;
        while (error_correct_flanking_kmers(r, ps) || error_correct_start(r, ps) || error_correct_end(r, ps)) {
            changes++;
            if (changes >= MAX_BASE_CORRECTIONS) {
                return changes;
            }
        }
        return changes;
    }
    private boolean error_correct_flanking_kmers(SAMRecord r, PackedSequence ps) {
        for (int i = 1; i < ps.length() - k; i++) {
            long leftKmer = ps.getKmer(i - 1, k);
            long leftTransform = collapseLookup.getOrDefault(leftKmer, SENTINEL_VALUE);
            if (leftTransform != SENTINEL_VALUE) {
                long rightKmer = ps.getKmer(i + 1, k);
                long rightTransform = collapseLookup.getOrDefault(rightKmer, SENTINEL_VALUE);
                if (rightTransform != SENTINEL_VALUE) {
                    // check they both agree on the base change being made
                    long basesMatching = KmerEncodingHelper.basesMatching(k - 2, leftTransform >> 4, rightTransform);
                    if (basesMatching == k - 2) {
                        ps.setKmer(leftTransform, i - 1, k);
                        return true;
                    }
                }
            }
        }
        return false;
    }
    private boolean error_correct_start(SAMRecord r, PackedSequence ps) {
        long kmer = ps.getKmer(0, k);
        long transform = collapseLookup.getOrDefault(kmer, SENTINEL_VALUE);
        if (transform != SENTINEL_VALUE) {
            // check we want to change the first or second base
            int basesDifferent = KmerEncodingHelper.basesDifference(k - 2, kmer, transform);
            boolean lowBytesSame = ((kmer & 15) != (transform & 15));
            if (basesDifferent == 0) {
                ps.setKmer(transform, 0, k);
                return true;
            }
        }
        return false;
    }
    private boolean error_correct_end(SAMRecord r, PackedSequence ps) {
        long kmer = ps.getKmer(ps.length() - k, k);
        long transform = collapseLookup.getOrDefault(kmer, SENTINEL_VALUE);
        if (transform != SENTINEL_VALUE) {
            // check we want to change the last or second last base
            if ((kmer & 15) != (transform & 15)) {
                ps.setKmer(transform, ps.length() - k, k);
                return true;
            }
        }
        return false;
    }

    private void ensureCollapseLookup() {
        if (collapseLookup == null) {
            collapseLookup = new Long2LongOpenHashMap();
            for (long kmer : kmerCounts.keySet()) {
                int count = kmerCounts.get(kmer);
                long bestNeighbourKmer = bestNeighbour(kmer);
                int bestNeighbourCount = kmerCounts.get(bestNeighbourKmer);
                if (count * collapseMultiple <= bestNeighbourCount) {
                    collapseLookup.put(kmer, bestNeighbourKmer);
                }
            }
            log.debug(String.format("Collapsed %d of %d kmer.", collapseLookup.size(), kmerCounts.size()));
        }
    }

    /**
     * Neighbour (hamming distance = 1) kmer with highest kmer count
     * @param kmer kmer to check neighbours of
     * @return kmer of highest neighbour. Returns the input kmer if all neighbours have 0 kmer counts.
     */
    private long bestNeighbour(long kmer) {
        long bestKmer = kmer;
        int bestCount = 0;
        for (long neighbour : KmerEncodingHelper.neighbouringStates(k, kmer)) {
            int count = kmerCounts.get(neighbour);
            if (count > bestCount) {
                bestKmer = (int)neighbour;
                bestCount = count;
            }
        }
        return bestKmer;
    }
}
