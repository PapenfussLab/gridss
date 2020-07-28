package au.edu.wehi.idsv.kraken;


import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

public class KrakenClassificationTest {
    @Test
    @Ignore("NYI")
    public void should_parse_lca() {
    }
    @Test
    @Ignore("NYI")
    public void should_parse_sequence_length() {
    }
    @Test
    public void should_parse_basic_fields() {
        KrakenClassification aligned = new KrakenClassification(new String("C\tid1\t1\t10\t1:10"));
        KrakenClassification unaligned = new KrakenClassification(new String("U\tid2\t0\t100\t"));
        Assert.assertTrue(aligned.isClassified);
        Assert.assertEquals("id1", aligned.sequenceId);
        Assert.assertEquals(1, aligned.taxonomyId);
        Assert.assertFalse(unaligned.isClassified);
        Assert.assertEquals("id2", unaligned.sequenceId);
        Assert.assertEquals(0, unaligned.taxonomyId);
    }
}