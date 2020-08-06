package gridss.kraken;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

public class AnnotateVariantsKrakenTest extends IntermediateFilesTest {
    @Test
    public void shouldAnnotate() throws IOException {
        VariantContext vc1 = new VariantContextBuilder()
                .start(100)
                .stop(100)
                .chr("polyA")
                .alleles("A", "AAAAAAAAAAAAAAAAAAAAATTTTTTTTTT.")
                .id("id1")
                .make();
        VariantContext vc2 = new VariantContextBuilder(vc1).id("id2").make();
        VariantContext vc3 = new VariantContextBuilder(vc1).id("id3").alleles("A", "TTTTT").make();
        VariantContext vc4 = new VariantContextBuilder(vc1).id("id4").make();
        VariantContext vc5 = new VariantContextBuilder(vc1).id("id5").alleles("A", "N[polyA:100[").make();
        VariantContext vc6 = new VariantContextBuilder(vc1).id("id6").alleles("A", "ATTTTTTTTTTTTTTTTTTTTTTTTT[polyA:100[").make();

        File vcf = new File(testFolder.getRoot(), "test.vcf");

        try (VariantContextWriter writer = new VariantContextWriterBuilder()
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .setOutputFile(vcf).build()) {
            VCFHeader header = new VCFHeader();
            header.setSequenceDictionary(SMALL_FA.getSequenceDictionary());
            writer.writeHeader(header);
            writer.add(vc1);
            writer.add(vc2);
            writer.add(vc3);
            writer.add(vc4);
            writer.add(vc5);
            writer.add(vc6);
        }
        File krakenOutput = new File(testFolder.getRoot(), "kraken.out");
        Files.write(krakenOutput.toPath(), ImmutableList.of(
                "C	id1	10	30	1:30",
                "U	id2	0	30	",
                "U	id4	0	30	",
                "C	id6	9606	25	9606:20"
        ), StandardCharsets.UTF_8);
        AnnotateVariantsKraken cmd = new AnnotateVariantsKraken();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = vcf;
        cmd.KRAKEN_INPUT = krakenOutput;
        cmd.OUTPUT = output;
        cmd.doWork();

        List<VariantContext> result = Lists.newArrayList(new VCFFileReader(output).iterator());

        Assert.assertEquals(10, result.get(0).getAttributeAsInt("INSTAXID", -1));
        Assert.assertFalse(result.get(1).hasAttribute("INSTAXID")); // No hit
        Assert.assertFalse(result.get(2).hasAttribute("INSTAXID")); // Not sent
        Assert.assertFalse(result.get(3).hasAttribute("INSTAXID")); // Not sent
        Assert.assertFalse(result.get(4).hasAttribute("INSTAXID")); // Not sent
        Assert.assertEquals(9606, result.get(5).getAttributeAsInt("INSTAXID", -1));
    }
}