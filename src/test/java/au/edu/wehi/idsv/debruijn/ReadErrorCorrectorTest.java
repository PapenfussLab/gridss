package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class ReadErrorCorrectorTest extends TestHelper {
    private static final String SEQ = S(RANDOM).substring(0, 100);
    private String intoSeq(int k, String seq, float threshold, int n) {
        List<SAMRecord> reads = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            SAMRecord r = withName("seq" + i, withSequence(B(SEQ), Read(0, 1, "100M")))[0];
            reads.add(r);
        }
        SAMRecord r = Read(0, 1, "100M");
        reads.add(r);
        r.setReadBases(B(seq));
        ReadErrorCorrector rec = new ReadErrorCorrector(k, threshold);
        reads.stream().forEach(x -> rec.countKmers(x));
        reads.stream().forEach(x -> rec.errorCorrect(x));
        return r.getReadString();
    }
    @Test
    public void should_correct_to_neighbour_sequence() {
        ReadErrorCorrector rec = new ReadErrorCorrector(21, 10);
        List<SAMRecord> reads = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            SAMRecord r = withSequence(B(SEQ), Read(0, 1, "100M"))[0];
            r.getReadBases()[i] = (byte)(r.getReadBases()[i] == 'T' ? 'A' : 'T');
            reads.add(r);
        }
        reads.stream().forEach(r -> rec.countKmers(r));
        reads.stream().forEach(r -> rec.errorCorrect(r));
        for (SAMRecord r : reads) {
            Assert.assertEquals(SEQ, S(r.getReadBases()));
        }
    }
    @Test
    public void should_respect_collapse_threshold() {
        String seq = "T" + SEQ.substring(1);
        Assert.assertEquals(seq, intoSeq(21, seq, 10, 9));
        Assert.assertNotEquals(seq, intoSeq(21, seq, 10, 10));
        Assert.assertNotEquals(seq, intoSeq(21, seq, 10, 11));
    }
    @Test
    public void should_correct_start() {
        String seq = "T" + SEQ.substring(1);
        Assert.assertEquals(SEQ, intoSeq(21, seq, 10, 20));
    }
    @Test
    public void should_correct_second() {
        String seq = "CT" + SEQ.substring(2);
        Assert.assertEquals(SEQ, intoSeq(21, seq, 10, 20));
    }
    @Test
    public void should_correct_third() {
        int offset = 2;
        String seq = SEQ.substring(0, offset) + "G" + SEQ.substring(offset + 1);
        Assert.assertEquals(SEQ, intoSeq(21, seq, 10, 20));
    }
    @Test
    public void should_correct_end() {
        String seq = SEQ.substring(0, SEQ.length() - 1) + "A";
        Assert.assertEquals(SEQ, intoSeq(21, seq, 10, 20));
    }
    @Test
    public void should_correct_second_last() {
        String seq = SEQ.substring(0, 98) + "TT";
        Assert.assertEquals(SEQ, intoSeq(21, seq, 10, 20));
    }
    @Test
    public void should_not_chain_extend_error_correction() {
        ReadErrorCorrector rec = new ReadErrorCorrector(4, 5);
        SAMRecord r = withSequence(B("AACTAAAAAAAAAAAAAAAAAAAACGTAACCGGTT"), Read(0, 1, "13M"))[0];
        rec.countKmers(r);
        rec.errorCorrect(r);
        Assert.assertEquals("AACTAAAAAAAAAAAAAAAAAAAACGTAACCGGTT", S(r.getReadBases()));
    }
    @Test
    public void should_not_correct_MNV() {
        StringBuilder sb = new StringBuilder(SEQ);
        sb.setCharAt(50, (char)SequenceUtil.complement((byte)sb.charAt(50)));
        sb.setCharAt(51, (char)SequenceUtil.complement((byte)sb.charAt(51)));
        sb.setCharAt(52, (char)SequenceUtil.complement((byte)sb.charAt(52)));
        Assert.assertEquals(SEQ, intoSeq(21, sb.toString(), 10, 20));
    }
}