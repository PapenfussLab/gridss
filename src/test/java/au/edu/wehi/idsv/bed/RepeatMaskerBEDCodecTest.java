package au.edu.wehi.idsv.bed;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.compress.utils.Lists;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

public class RepeatMaskerBEDCodecTest extends TestHelper {
    @Test
    public void should_parse_as_RMBedFeatures() throws IOException {
        List<BEDFeature> records = Lists.newArrayList(AbstractFeatureReader.getFeatureReader(
                new File("src/test/resources/hg19.rm.bedops.bed").getPath(),
                new RepeatMaskerBEDCodec(), false).iterator());
        Assert.assertTrue(records.stream().allMatch(x -> x instanceof RepeatMaskerBEDFeature));
    }

    /**
     * https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html
     */
    @Test
    public void should_parse_according_to_bedops_rmsk2bed_definition() throws IOException {
        List<RepeatMaskerBEDFeature> records = Lists.newArrayList(AbstractFeatureReader.getFeatureReader(
                new File("src/test/resources/hg19.rm.bedops.bed").getPath(),
                new RepeatMaskerBEDCodec(), false).iterator())
                .stream()
                .map(x -> (RepeatMaskerBEDFeature) x).collect(Collectors.toList());
        Assert.assertEquals(10, records.size());
        Assert.assertEquals(10001, records.get(0).getStart());
        Assert.assertEquals(10468, records.get(0).getEnd());
        Assert.assertEquals("(CCCTAA)n", records.get(0).getRepeatType());
        Assert.assertEquals("Simple_repeat", records.get(0).getRepeatClass());
        Assert.assertEquals(Strand.POSITIVE, records.get(0).getStrand());
        Assert.assertEquals(10469, records.get(1).getStart());
        Assert.assertEquals(11447, records.get(1).getEnd());
        Assert.assertEquals("TAR1", records.get(1).getRepeatType());
        Assert.assertEquals("Satellite/telo", records.get(1).getRepeatClass());
        Assert.assertEquals(Strand.NEGATIVE, records.get(1).getStrand());
    }
    @Test
    public void should_cache_strings() throws IOException {
        List<RepeatMaskerBEDFeature> records = Lists.newArrayList(AbstractFeatureReader.getFeatureReader(
                new File("src/test/resources/hg19.rm.bedops.bed").getPath(),
                new RepeatMaskerBEDCodec(), false).iterator())
                .stream()
                .map(x -> (RepeatMaskerBEDFeature) x).collect(Collectors.toList());
        Assert.assertTrue(records.get(0).getContig() == records.get(1).getContig());
    }
}