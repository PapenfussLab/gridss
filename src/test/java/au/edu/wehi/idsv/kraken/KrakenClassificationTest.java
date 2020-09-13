package au.edu.wehi.idsv.kraken;


import org.junit.Assert;
import org.junit.Test;

public class KrakenClassificationTest {
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
    @Test
    public void should_parse_length() {
        KrakenClassification aligned1 = new KrakenClassification(new String("C\tid2\t9606\t15\t"));
        KrakenClassification aligned2 = new KrakenClassification(new String("C\tid1\t1\t10|20\t1:10"));
        Assert.assertEquals(15, aligned1.sequenceLength);
        Assert.assertEquals(0, aligned1.sequenceLength2);
        Assert.assertEquals(10, aligned2.sequenceLength);
        Assert.assertEquals(20, aligned2.sequenceLength2);
    }
    @Test
    public void should_parse_kmer_alignments_se() {
        KrakenClassification kc = new KrakenClassification(new String("C\tid1\t1\t10|20\t1:10"));
        Assert.assertEquals(1, kc.kmerTaxonomyIds.size());
        Assert.assertEquals(10, kc.kmerTaxonomyIds.get(0).kmerCount);
        Assert.assertEquals(1, kc.kmerTaxonomyIds.get(0).taxonomyId);
    }
    @Test
    public void should_parse_kmer_alignments_pe() {
        KrakenClassification kc = new KrakenClassification(new String("C\tid2\t9606\t15\t562:13 561:4 A:31 0:1 562:3|:|9606:15"));
        Assert.assertEquals(5, kc.kmerTaxonomyIds.size());
        Assert.assertEquals(1, kc.kmerTaxonomyIds2.size());
        Assert.assertEquals(13, kc.kmerTaxonomyIds.get(0).kmerCount);
        Assert.assertEquals(4, kc.kmerTaxonomyIds.get(1).kmerCount);
        Assert.assertEquals(31, kc.kmerTaxonomyIds.get(2).kmerCount);
        Assert.assertEquals(1, kc.kmerTaxonomyIds.get(3).kmerCount);
        Assert.assertEquals(3, kc.kmerTaxonomyIds.get(4).kmerCount);
        Assert.assertEquals(15, kc.kmerTaxonomyIds2.get(0).kmerCount);
        Assert.assertEquals(562, kc.kmerTaxonomyIds.get(0).taxonomyId);
        Assert.assertEquals(561, kc.kmerTaxonomyIds.get(1).taxonomyId);
        Assert.assertEquals(KrakenKmerClassification.AMBIGUOUS, kc.kmerTaxonomyIds.get(2).taxonomyId);
        Assert.assertEquals(0, kc.kmerTaxonomyIds.get(3).taxonomyId);
        Assert.assertEquals(562, kc.kmerTaxonomyIds.get(4).taxonomyId);
        Assert.assertEquals(9606, kc.kmerTaxonomyIds2.get(0).taxonomyId);
    }
}