package au.edu.wehi.idsv.sim;

import java.util.Arrays;
import java.util.Random;

public class RandomReferenceSequenceGenerator implements RandomSequenceGenerator {
    /**
     * Maximum number of subsequences to sample before giving up because they all have an N.
     */
    private static final int MAXIMUM_ATTEMPTS = 128;
    private static final boolean ALLOW_N = false;
    private final Random random;
    private final byte[] seq;
    private int lastStartPos;
    private int lastEndPos;

    public RandomReferenceSequenceGenerator(int seed, byte[] seq) {
        this.random = new Random(seed);
        this.seq = seq;
    }
    public byte[] getBases(int length) {
        if (length > seq.length) throw new IllegalArgumentException("sequence cannot be longer than reference contig length when sampling from the reference contig");
        byte[] subseq = null;
        for (int i = 0; i < MAXIMUM_ATTEMPTS; i++) {
            int start = random.nextInt(seq.length - length + 1);
            subseq = Arrays.copyOfRange(seq, start, start + length);
            if (!ALLOW_N) {
                for (int j = 0; j < subseq.length; j++) {
                    if (subseq[j] == 'N') {
                        // No Ns allowed
                        continue;
                    }
                }
            }
            this.lastStartPos = start;
            this.lastEndPos = start + subseq.length - 1;
            break;
        }
        return subseq;
    }
    public Integer lastStartPosition() { return lastStartPos; }
    public Integer lastEndPosition() { return lastEndPos; }
}
