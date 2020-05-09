package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

public class VariantContextRepeatMaskerAnnotatorTest extends TestHelper {
    @Test
    public void should_match_with_repeat_masker_interval() throws IOException {
        VariantContextRepeatMaskerAnnotator rma = new VariantContextRepeatMaskerAnnotator(new File("src/test/resources/hg19.rm.bedops.bed"));
        VariantContext vc = minimalVariant().attribute("BEALN", "chr1:9999|-|4M|60").make();
        vc = rma.apply(vc);
        Assert.assertEquals("(CCCTAA)n", vc.getAttributeAsString("INSRMRT", null));
        Assert.assertEquals("Simple_repeat", vc.getAttributeAsString("INSRMRC", null));
        Assert.assertEquals("-", vc.getAttributeAsString("INSRMRO", null));
        Assert.assertEquals(0.5, vc.getAttributeAsDouble("INSRMP", 0), 0);
    }
    @Test
    public void should_match_overlap_bounds_as_closed_interval() throws IOException {
        VariantContextRepeatMaskerAnnotator rma = new VariantContextRepeatMaskerAnnotator(new File("src/test/resources/hg19.rm.bedops.bed"));
        Assert.assertFalse(rma.apply(minimalVariant().attribute("BEALN", "chr1:9999|-|1M|60").make()).hasAttribute("INSRMRT"));
        Assert.assertFalse(rma.apply(minimalVariant().attribute("BEALN", "chr1:10000|-|1M|60").make()).hasAttribute("INSRMRT"));
        Assert.assertTrue(rma.apply(minimalVariant().attribute("BEALN", "chr1:10001|-|1M|60").make()).hasAttribute("INSRMRT"));
        Assert.assertTrue(rma.apply(minimalVariant().attribute("BEALN", "chr1:22074|-|1M|60").make()).hasAttribute("INSRMRT"));
        Assert.assertTrue(rma.apply(minimalVariant().attribute("BEALN", "chr1:22075|-|1M|60").make()).hasAttribute("INSRMRT"));
        Assert.assertFalse(rma.apply(minimalVariant().attribute("BEALN", "chr1:22076|-|1M|60").make()).hasAttribute("INSRMRT"));
    }
    @Test
    public void should_match_with_greatest_overlap_for_alignment() throws IOException {
        VariantContextRepeatMaskerAnnotator rma = new VariantContextRepeatMaskerAnnotator(new File("src/test/resources/hg19.rm.bedops.bed"));
        VariantContext vc = minimalVariant().attribute("BEALN", ImmutableList.of("chr1:10465|-|10M|60")).make();
        vc = rma.apply(vc);
        Assert.assertEquals("TAR1", vc.getAttributeAsString("INSRMRT", null));
    }
    @Test
    public void should_match_with_greatest_overlap_for_any_alignment() throws IOException {
        VariantContextRepeatMaskerAnnotator rma = new VariantContextRepeatMaskerAnnotator(new File("src/test/resources/hg19.rm.bedops.bed"));
        VariantContext vc = minimalVariant().attribute("BEALN", ImmutableList.of("chr1:10465|-|10M|60", "chr1:11678|-|100M|60")).make();
        vc = rma.apply(vc);
        Assert.assertEquals("MER5B", vc.getAttributeAsString("INSRMRT", null));
        Assert.assertEquals("+", vc.getAttributeAsString("INSRMRO", null));
    }
}