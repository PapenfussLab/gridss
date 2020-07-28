package au.edu.wehi.idsv.kraken;

import org.apache.commons.compress.utils.Lists;
import org.junit.Assert;
import org.junit.Test;

import java.io.*;
import java.util.List;

public class KrakenParserTest {
    @Test
    public void should_parse_kraken2_output() throws FileNotFoundException {
        List<KrakenClassification> result = Lists.newArrayList(new KrakenParser(new BufferedReader(new InputStreamReader(new FileInputStream(new File("src/test/resources/kraken2_output.tsv"))))));
        Assert.assertEquals(5, result.size());

        Assert.assertFalse(result.get(0).isClassified);
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1245:9625:11303", result.get(0).sequenceId);
        Assert.assertEquals(0, result.get(0).taxonomyId);

        Assert.assertFalse(result.get(1).isClassified);
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1245:9625:11303", result.get(1).sequenceId);
        Assert.assertEquals(0, result.get(1).taxonomyId);

        Assert.assertTrue(result.get(2).isClassified);
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1245:9670:18239", result.get(2).sequenceId);
        Assert.assertEquals(9606, result.get(2).taxonomyId);

        Assert.assertTrue(result.get(3).isClassified);
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1245:9670:18239", result.get(3).sequenceId);
        Assert.assertEquals(9606, result.get(3).taxonomyId);

        Assert.assertTrue(result.get(4).isClassified);
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2146:29749:32377", result.get(4).sequenceId);
        Assert.assertEquals(28384, result.get(4).taxonomyId);
    }
}