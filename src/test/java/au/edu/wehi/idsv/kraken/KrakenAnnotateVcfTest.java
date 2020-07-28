package au.edu.wehi.idsv.kraken;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.junit.Assert;
import org.junit.Test;
import samtools.htsjdk.fastq.MockFastqWriter;

import java.io.*;
import java.util.ArrayList;

public class KrakenAnnotateVcfTest {
    @Test
    public void should_annotate_vcf() {
        VariantContext vc1 = new VariantContextBuilder()
                .start(100)
                .stop(100)
                .chr("chr1")
                .alleles("A", "AAAAAAAAAAAAAAAAAAAAATTTTTTTTTT.")
                .id("id1")
                .make();
        VariantContext vc2 = new VariantContextBuilder(vc1).id("id2").make();
        VariantContext vc3 = new VariantContextBuilder(vc1).id("id3").alleles("A", "TTTTT").make();
        VariantContext vc4 = new VariantContextBuilder(vc1).id("id4").make();
        VariantContext vc5 = new VariantContextBuilder(vc1).id("id5").alleles("A", "N[chr2:100[").make();
        VariantContext vc6 = new VariantContextBuilder(vc1).id("id6").alleles("A", "ATTTTTTTTTTTTTTTTTTTTTTTTT[chr2:100[").make();

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(baos);
        ps.println("C	id1	10	30	1:30");
        ps.println("U	id2	0	30	");
        ps.println("U	id4	0	30	");
        ps.println("C	id6	9606	25	9606:20");
        ps.close();

        KrakenParser parser = new KrakenParser(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(baos.toByteArray()))));
        MockFastqWriter fqw = new MockFastqWriter();
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", 1000));
        dict.addSequence(new SAMSequenceRecord("chr2", 1000));
        KrakenAnnotateVcf kav = new KrakenAnnotateVcf(parser, fqw, ImmutableList.of(vc1, vc2, vc3, vc4, vc5, vc6).iterator(), dict, 10);


        ArrayList<VariantContext> result = Lists.newArrayList(kav);
        Assert.assertEquals(4, fqw.records.size());
        Assert.assertEquals(10, result.get(0).getAttribute("INSTAXID"));
        Assert.assertFalse(result.get(1).hasAttribute("INSTAXID")); // No hit
        Assert.assertFalse(result.get(2).hasAttribute("INSTAXID")); // Not sent
        Assert.assertFalse(result.get(3).hasAttribute("INSTAXID")); // Not sent
        Assert.assertFalse(result.get(4).hasAttribute("INSTAXID")); // Not sent
        Assert.assertEquals(9606, result.get(5).getAttribute("INSTAXID"));

    }
}