package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.configuration.ErrorCorrectionConfiguration;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import it.unimi.dsi.fastutil.Hash;
import it.unimi.dsi.fastutil.longs.Long2IntMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import it.unimi.dsi.fastutil.longs.LongSet;

import java.util.HashSet;
import java.util.Set;

public class ReadErrorCorrector {
    private static final Log log = Log.getInstance(ReadErrorCorrector.class);
    private final Long2IntMap kmerCounts = new Long2IntOpenHashMap();
    private final int k;
    private final float kmerErrorCorrectionMultiple;
    private final boolean deduplicateReadKmers;
    private final int maxCorrectionsInKmer;
    private int maxCount = 0;
    private int maxCollapseCount = 0;

    public ReadErrorCorrector(ErrorCorrectionConfiguration ecc) {
        this(ecc.k, ecc.kmerErrorCorrectionMultiple, ecc.maxCorrectionsInKmer, ecc.deduplicateReadKmers);
    }
    public ReadErrorCorrector(int k, float kmerErrorCorrectionMultiple, int maxCorrectionsInKmer, boolean deduplicateReadKmers) {
        if (k > 31) throw new IllegalArgumentException("k cannot exceed 31");
        this.k = k;
        this.kmerErrorCorrectionMultiple = kmerErrorCorrectionMultiple;
        this.maxCorrectionsInKmer = maxCorrectionsInKmer;
        this.deduplicateReadKmers = deduplicateReadKmers;
    }

    /**
     * Error corrects the reads underlying the given evidence.
     * @param k kmer size
     * @param kmerErrorCorrectionMultiple kmer abundance relative to maximum kmer before collapsing. For example, if the maximum kmer abundance is 100, kmers found less than 10 times are error correction candidates.
     * @param deduplicateReadKmers determines whether read kmers are counted based on the abundance or just presence/absence
     *                             Counting based on abundance causes issues when homopolymers runs are present.
     * @param evidence evidence to perform error correction on
     */
    public static void errorCorrect(int k, float kmerErrorCorrectionMultiple, int maxCorrectionsInKmer, boolean deduplicateReadKmers, Iterable<? extends DirectedEvidence> evidence) {
        ReadErrorCorrector ec = new ReadErrorCorrector(k, kmerErrorCorrectionMultiple, maxCorrectionsInKmer, deduplicateReadKmers);
        // need to deduplicate the underlying reads so we don't double count
        // kmers from reads with multiple evidence (e.g. multiple indels or SC on both ends)
        Set<SAMRecord> reads = new HashSet<>();
        Set<SAMRecord> rcreads = new HashSet<>();
        for (DirectedEvidence de : evidence) {
            reads.add(de.getUnderlyingSAMRecord());
            if (de instanceof NonReferenceReadPair) {
                SAMRecord mate = ((NonReferenceReadPair) de).getNonReferenceRead();
                if ((de.getBreakendSummary().direction == BreakendDirection.Forward) ^ mate.getReadNegativeStrandFlag()) {
                    rcreads.add(mate);
                } else {
                    reads.add(mate);
                }
            }
        }
        // TODO: cache PackedSequence
        reads.stream().forEach(r -> ec.countKmers(r, false));
        rcreads.stream().forEach(r -> ec.countKmers(r, true));
        reads.stream().forEach(r -> ec.errorCorrect(r, false));
        rcreads.stream().forEach(r -> ec.errorCorrect(r, true));
    }

    public void countKmers(SAMRecord r, boolean reverseComplement) {
        PackedSequence ps = new PackedSequence(r.getReadBases(), reverseComplement, reverseComplement);
        LongSet encountered = new LongOpenHashSet(Math.max(1, ps.length() - k + 1), Hash.FAST_LOAD_FACTOR);
        for (int i = 0; i < ps.length() - k + 1; i++) {
            long kmer = ps.getKmer(i, k);
            if (this.deduplicateReadKmers && encountered.contains(kmer)) {
                continue;
            }
            int count = kmerCounts.get(kmer);
            count++;
            kmerCounts.put(kmer, count);
            if (count > maxCount) {
                maxCount = count;
                refreshMaxCollapseCount();
            }
            if (this.deduplicateReadKmers) {
                encountered.add(kmer);
            }
        }
    }
    public int errorCorrect(SAMRecord r, boolean reverseComplement) {
        if (r.getReadLength() < k) return 0;
        PackedSequence ps = new PackedSequence(r.getReadBases(), reverseComplement, reverseComplement);
        PackedSequence oldps = new PackedSequence(ps);
        int basesChanged = 0;
        basesChanged += musket_two_sided(ps);
        basesChanged += musket_one_sided_greedy_without_voting(ps);
        if (basesChanged > 0) {
            if (basesChanged > maxCorrectionsInKmer && tooManyDifferencesInWindow(ps, oldps)) {
                // discard the changes - too much in a window
                basesChanged = 0;
            } else {
                byte[] seq = ps.getBytes(0, r.getReadBases().length);
                if (reverseComplement) {
                    SequenceUtil.reverseComplement(seq);
                }
                r.setReadBases(seq);
            }
        }
        //debug_dump_changes(r, ps, oldps);
        return basesChanged;
    }

    private boolean tooManyDifferencesInWindow(PackedSequence ps1, PackedSequence ps2) {
        for (int i = 0; i < ps1.length() - (k - 1); i++) {
            int diff = KmerEncodingHelper.basesDifference(k, ps1.getKmer(i, k), ps2.getKmer(i, k));
            if (diff > maxCorrectionsInKmer) {
                return true;
            }
        }
        return false;
    }

    private void debug_dump_changes(SAMRecord r, PackedSequence ps, PackedSequence oldps) {
        String oldSeq = new String(oldps.getBytes(0, r.getReadBases().length));
        String newSeq = new String(ps.getBytes(0, r.getReadBases().length));
        System.err.printf("\n%s\n", oldSeq);
        for (int i = 0; i < oldSeq.length() - (k - 1); i++) {
            if (isSafeBase(oldps, i)) {
                System.err.printf("*");
            } else {
                System.err.printf(" ");
            }
        }
        System.err.printf("\n");
        for (int i = 0; i < oldSeq.length(); i++) {
            if (oldSeq.charAt(i) == newSeq.charAt(i)) {
                System.err.printf(" ");
            } else {
                System.err.printf("|");
            }
        }
        System.err.printf("\n");
        for (int i = 0; i < oldSeq.length(); i++) {
            if (oldSeq.charAt(i) == newSeq.charAt(i)) {
                System.err.printf(" ");
            } else {
                System.err.printf("%s", newSeq.charAt(i));
            }
        }
        System.err.printf("\n");
    }

    /**
     * Max count that we will consider collapsing. All kmers above this threshold are safe from
     * error correction
     */
    private void refreshMaxCollapseCount() {
        maxCollapseCount = (int)Math.floor(maxCount / kmerErrorCorrectionMultiple);
    }

    private int musket_two_sided(PackedSequence ps) {
        //List<Integer> zzcounts = IntStream.range(0, ps.length() - k + 1).mapToObj(offset -> kmerCounts.get(ps.getKmer(offset, k))).collect(Collectors.toList());
        int i = k - 1;
        int basesChanged = 0;
        while (i + k < ps.length()) {
            long rightKmer = ps.getKmer(i, k);
            int rightCount = kmerCounts.get(rightKmer);
            if (rightCount > maxCollapseCount) {
                i += k; // advance to first non-overlapping
                continue;
            }
            int leftKmerOffset = i - (k - 1);
            long leftKmer = ps.getKmer(leftKmerOffset, k);
            int leftCount = kmerCounts.get(leftKmer);
            if (leftCount > maxCollapseCount) {
                i += 1; // advance to first non-overlapping kmer
                continue;
            }
            // check both kmers for their best replacement
            long leftCollapse = neighbourToCollapseInto(leftKmer, k - 1, leftCount);
            if (leftCollapse == leftKmer) {
                i++;
                continue;
            }
            long rightCollapse = neighbourToCollapseInto(rightKmer, 0, rightCount);
            if (rightCollapse == rightKmer) {
                i++;
                continue;
            }
            if (kmerCounts.get(leftCollapse) >= kmerCounts.get(rightCollapse)) {
                // collapse using replacement base from left kmer
                ps.setKmer(leftCollapse, leftKmerOffset, k);
            } else {
                ps.setKmer(rightCollapse, i, k);
            }
            i += k; // right kmer is now good
            basesChanged++;
        }
        return basesChanged;
    }
    private int musket_one_sided_greedy_without_voting(PackedSequence ps) {
        int safePosition = advanceToSafePosition(ps, 0);
        int basesChanged = 0;
        // nothing is safe
        if (safePosition > ps.length() - k) return basesChanged;
        if (safePosition > 0) {
            basesChanged += fixBaseAdjacentTo(ps, safePosition, -1);
        }
        safePosition = advanceToEndOfSafeRegion(ps, safePosition);
        while (safePosition <= ps.length() - k - 1) {
            basesChanged += fixBaseAdjacentTo(ps, safePosition, 1);
            safePosition = advanceToSafePosition(ps, safePosition + 1);
            safePosition = advanceToEndOfSafeRegion(ps, safePosition);
        }
        return basesChanged;
    }

    private int advanceToSafePosition(PackedSequence ps, int offset) {
        while (offset <= ps.length() - k) {
            if (isSafeKmerStartPosition(ps, offset)) {
                return offset;
            }
            offset++;
        }
        return offset;
    }
    private boolean isSafeKmerStartPosition(PackedSequence ps, int offset) {
        return kmerCounts.get(ps.getKmer(offset, k)) > maxCollapseCount;
    }
    private boolean isSafeBase(PackedSequence ps, int offset) {
        for (int i = Math.max(0, offset - (k - 1)); i <= Math.min(offset, ps.length() - k); i++) {
            if (isSafeKmerStartPosition(ps, i)) {
                return true;
            }
        }
        return false;
    }

    private int advanceToEndOfSafeRegion(PackedSequence ps, int offset) {
        while (offset <= ps.length() - k) {
            if (!isSafeKmerStartPosition(ps, offset)) {
                // unsafe kmer - we've hit the end of our safe region
                return offset - 1;
            }
            int strideOffset = offset + k;
            if (strideOffset <= ps.length() - k && isSafeKmerStartPosition(ps, strideOffset)) {
                // we can jump forward if the next non-overlapping kmer is also safe
                offset = strideOffset;
            } else {
                offset++;
            }
        }
        return offset;
    }

    private int fixBaseAdjacentTo(PackedSequence ps, int position, int direction) {
        int nextPosition = position + direction;
        // out of bounds
        if (nextPosition < 0 || nextPosition > ps.length() - k) return 0;
        int nextNextPosition = nextPosition + direction;
        long nextKmer = ps.getKmer(nextPosition, k);
        int nextKmerBaseOfInterest = direction == -1 ? 0 : k - 1;
        long nextCollapseKmer = neighbourToCollapseInto(nextKmer, nextKmerBaseOfInterest, kmerCounts.get(nextKmer));
        if (nextCollapseKmer == nextKmer) return 0;
        if (nextNextPosition < 0 || nextNextPosition > ps.length() - k) {
            // we're at the start or end of the read so we can't check the next base
            ps.setKmer(nextCollapseKmer, nextPosition, k);
            return 1;
        }
        // check that the next kmer also agrees with the error correction
        long nextNextKmer = ps.getKmer(nextNextPosition, k);
        int nextNextKmerBaseOfInterest = direction == -1 ? 1 : k - 2;
        long nextNextCollapseKmer = neighbourToCollapseInto(nextNextKmer, nextNextKmerBaseOfInterest, kmerCounts.get(nextNextKmer));

        long nextKmerChangedBase = KmerEncodingHelper.getBase(k, nextCollapseKmer, nextKmerBaseOfInterest);
        long nextNextChangedBase = KmerEncodingHelper.getBase(k, nextNextCollapseKmer,nextNextKmerBaseOfInterest);
        // make sure they agree on the base to be changed
        if (nextKmerChangedBase != nextNextChangedBase) return 0;
        ps.setKmer(nextCollapseKmer, nextPosition, k);
        return 1;
    }

    private long neighbourToCollapseInto(long kmer, int baseOffset, int count) {
        long bestKmer = kmer;
        // use 1 less than the collapse threshold as the sentinel value
        int bestCount = (int)Math.ceil(count * kmerErrorCorrectionMultiple) - 1;
        for (long j = 1; j < 4; j++) { // XOR 0 = self so we can start at 1
            long neighbourKmer = kmer ^ (j << ((k - 1 - baseOffset) * 2));
            int neighbourCount = kmerCounts.get(neighbourKmer);
            if (neighbourCount > bestCount) {
                bestKmer = neighbourKmer;
                bestCount = neighbourCount;
            }
        }
        return bestKmer;
    }
}
