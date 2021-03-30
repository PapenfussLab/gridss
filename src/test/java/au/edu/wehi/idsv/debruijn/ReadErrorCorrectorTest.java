package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.DirectedEvidence;
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
        reads.stream().forEach(x -> rec.countKmers(x, false));
        rec.errorCorrect(r, false);
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
        reads.stream().forEach(r -> rec.countKmers(r, false));
        reads.stream().forEach(r -> rec.errorCorrect(r, false));
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
    public void should_correct_middle() {
        StringBuilder sb = new StringBuilder(SEQ);
        sb.setCharAt(50, (char)SequenceUtil.complement((byte)sb.charAt(50)));
        String seq = sb.toString();
        Assert.assertEquals(SEQ, intoSeq(21, seq, 10, 20));
    }
    @Test
    public void should_not_chain_extend_error_correction() {
        ReadErrorCorrector rec = new ReadErrorCorrector(4, 5);
        SAMRecord r = withSequence(B("AACTAAAAAAAAAAAAAAAAAAAACGTAACCGGTT"), Read(0, 1, "13M"))[0];
        rec.countKmers(r, false);
        rec.errorCorrect(r, false);
        Assert.assertEquals("AACTAAAAAAAAAAAAAAAAAAAACGTAACCGGTT", S(r.getReadBases()));
    }
    @Test
    public void should_not_correct_MNV() {
        StringBuilder sb = new StringBuilder(SEQ);
        sb.setCharAt(50, (char)SequenceUtil.complement((byte)sb.charAt(50)));
        sb.setCharAt(51, (char)SequenceUtil.complement((byte)sb.charAt(51)));
        sb.setCharAt(52, (char)SequenceUtil.complement((byte)sb.charAt(52)));
        String s2 = sb.toString();
        Assert.assertEquals(s2, intoSeq(21, s2, 10, 20));
    }
    @Test
    public void should_correct_mate_evidence_according_to_expected_orientation() {
        MockSAMEvidenceSource ses = SES();
        List<DirectedEvidence> evidence = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            SAMRecord r = withName("seq" + i, withSequence(B(SEQ), Read(2, 1, "50M50S")))[0];
            evidence.add(SCE(FWD, ses, r));
        }
        SAMRecord[] dpfp = DP(0, 2, "100M", true, 2, 1, "100M", true);
        dpfp[0].setReadBases(B("G" + SEQ.substring(1)));
        dpfp[1].setReadBases(B(SequenceUtil.reverseComplement("T" + SEQ.substring(1))));
        evidence.add(NRRP(ses, dpfp));

        SAMRecord[] dpfn = DP(0, 3, "100M", true, 2, 1, "100M", false);
        dpfn[0].setReadBases(B("G" + SEQ.substring(1)));
        dpfn[1].setReadBases(B("T" + SEQ.substring(1)));
        evidence.add(NRRP(ses, dpfn));

        SAMRecord[] dprp = DP(0, 4, "100M", false, 2, 1, "100M", false);
        dprp[0].setReadBases(B("G" + SEQ.substring(1)));
        dprp[1].setReadBases(B(SequenceUtil.reverseComplement("T" + SEQ.substring(1))));
        evidence.add(NRRP(ses, dprp));

        SAMRecord[] dprn = DP(0, 5, "100M", false, 2, 1, "100M", true);
        dprn[0].setReadBases(B("G" + SEQ.substring(1)));
        dprn[1].setReadBases(B("T" + SEQ.substring(1)));
        evidence.add(NRRP(ses, dprn));

        ReadErrorCorrector.errorCorrect(21, 5, evidence);
        Assert.assertEquals(SEQ, dpfp[0].getReadString());
        Assert.assertEquals(SequenceUtil.reverseComplement(SEQ), dpfp[1].getReadString());
        Assert.assertEquals(SEQ, dpfn[0].getReadString());
        Assert.assertEquals(SEQ, dpfn[1].getReadString());
        Assert.assertEquals(SEQ, dprp[0].getReadString());
        Assert.assertEquals(SequenceUtil.reverseComplement(SEQ), dprp[1].getReadString());
        Assert.assertEquals(SEQ, dprn[0].getReadString());
        Assert.assertEquals(SEQ, dprn[1].getReadString());
    }
    @Test
    public void should_handle_any_read_length() {
        MockSAMEvidenceSource ses = SES();
        List<DirectedEvidence> evidence = new ArrayList<>();
        for (int i = 1; i < 150; i++) {
            evidence.add(SCE(FWD, ses, Read(2, 1, "1M" + i + "S")));
            for (int j = 1; j < 150; j++) {
                evidence.add(NRRP(ses, DP(0, 1, i + "M", false, 1, 1, j + "M", true)));
            }
        }
        ReadErrorCorrector.errorCorrect(21, 5, evidence);
    }
}
