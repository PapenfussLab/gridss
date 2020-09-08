package gridss.repeatmasker;

import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class AnnotateVariantsRepeatMaskerTest extends IntermediateFilesTest {
    @Test
    public void should_annotate_rm_output() throws IOException {
        AnnotateVariantsRepeatMasker cmd = new AnnotateVariantsRepeatMasker();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/repeatmasker/ebv/merged.bwa.bam.gridss.vcf");
        cmd.REPEAT_MASKER = new File("src/test/resources/repeatmasker/ebv/merged.bwa.fa.out");
        cmd.OUTPUT = output;
        cmd.doWork();

        List<VariantContext> vcf = getRawVcf(output);

        VariantContext vc = vcf.stream().filter(x -> x.getID().equals("gridss312b_2342b")).findAny().get();
        Assert.assertEquals(79 / 80.0, vc.getAttributeAsDouble("INSRMP", 0), 0.001);
        Assert.assertEquals("rRNA", vc.getAttributeAsString("INSRMRC", ""));
        Assert.assertEquals("+", vc.getAttributeAsString("INSRMRO", ""));
        Assert.assertEquals("LSU-rRNA_Hsa", vc.getAttributeAsString("INSRMRT", ""));
        Assert.assertEquals("LSU-rRNA_Hsa#rRNA|2099|+|79M||", vc.getAttributeAsStringList("INSRM", "").get(0));
    }
    @Test
    public void should_annotate_rm_alignment() throws IOException {
        AnnotateVariantsRepeatMasker cmd = new AnnotateVariantsRepeatMasker();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/repeatmasker/ebv/merged.bwa.bam.gridss.vcf");
        cmd.REPEAT_MASKER = new File("src/test/resources/repeatmasker/ebv/merged.bwa.fa.align");
        cmd.OUTPUT = output;
        cmd.doWork();

        List<VariantContext> vcf = getRawVcf(output);

        VariantContext vc = vcf.stream().filter(x -> x.getID().equals("gridss312b_2342b")).findAny().get();
        Assert.assertEquals(79 / 80.0, vc.getAttributeAsDouble("INSRMP", 0), 0.001);
        Assert.assertEquals("rRNA", vc.getAttributeAsString("INSRMRC", ""));
        Assert.assertEquals("+", vc.getAttributeAsString("INSRMRO", ""));
        Assert.assertEquals("LSU-rRNA_Hsa", vc.getAttributeAsString("INSRMRT", ""));
        Assert.assertEquals("LSU-rRNA_Hsa#rRNA|2099|+|7=2M1=1I3=2X1M42=1M1=2I16=||5", vc.getAttributeAsStringList("INSRM", "").get(0));
    }
}