package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class VirusBreakendFilterTest extends IntermediateFilesTest {
    @Test
    public void should_filter_no_host_alignment() throws IOException {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("no_host_alignment")));
    }
    @Test
    public void should_filter_repeat_50_percent_overlap() throws IOException {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("repeat_overlap")));
    }
    @Test
    public void should_filter_short_host_hits() throws IOException {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("short_host_alignment")));
    }
    @Test
    public void should_filter_breakpoints() throws IOException {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("breakpoint")));
    }
    @Test
    public void should_retain_small_contained_repeat() throws IOException {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertTrue(vcf.stream().anyMatch(vc -> vc.getID().equals("small_contained_repeat")));
    }
    @Test
    public void should_retain_simple_repeat_no_overlap() throws IOException {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertTrue(vcf.stream().anyMatch(vc -> vc.getID().equals("repeat_no_overlap")));
    }
    @Test
    public void should_retain_complex_repeat_overlap() throws IOException {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertTrue(vcf.stream().anyMatch(vc -> vc.getID().equals("line_overlap")));
    }
}