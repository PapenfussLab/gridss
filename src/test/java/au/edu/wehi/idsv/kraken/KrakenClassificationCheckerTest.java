package au.edu.wehi.idsv.kraken;


import com.google.common.collect.ImmutableList;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

public class KrakenClassificationCheckerTest {
    @Test
    public void shouldIncludeSplitReads() throws IOException {
        KrakenClassificationChecker kkc = new KrakenClassificationChecker(ImmutableList.of(207598), new File("src/test/resources/ncbi/homo_sapiens.nodes.dmp"));

        // all ancestor not fine
        Assert.assertFalse(kkc.isOfInterest(new KrakenClassification("C\tid1\t1\t10\t1:1")));
        // good then by bad is fine
        Assert.assertTrue(kkc.isOfInterest(new KrakenClassification("C\tid1\t1\t10\t9606:1 10239:1")));
        // bad then good is fine
        Assert.assertTrue(kkc.isOfInterest(new KrakenClassification("C\tid1\t1\t10\t10239:1 9606:1")));
        // good flanked by bad is not fine
        Assert.assertFalse(kkc.isOfInterest(new KrakenClassification("C\tid1\t1\t10\t10239:1 9606:1 10239:1")));
        // unless kraken call it as good
        Assert.assertTrue(kkc.isOfInterest(new KrakenClassification("C\tid1\t9606\t10\t10239:1 9606:1 10239:1")));
        // good flanked by amb and ancestral is fine
        Assert.assertTrue(kkc.isOfInterest(new KrakenClassification("C\tid1\t10239\t10\tA:1 1:1 2759:1 9606:1 10239:1")));
        Assert.assertTrue(kkc.isOfInterest(new KrakenClassification("C\tid1\t10239\t10\t10239:1 9606:1 A:1 1:1 2759:1")));
    }

    @Test
    public void should_check_both_reads_in_split_read() throws IOException {
        KrakenClassificationChecker kkc = new KrakenClassificationChecker(ImmutableList.of(207598), new File("src/test/resources/ncbi/homo_sapiens.nodes.dmp"));
        // read2 is good
        Assert.assertTrue(kkc.isOfInterest(new KrakenClassification("C\tid1\t1\t10\t1:1|:|9606:1 10239:1")));
    }
}