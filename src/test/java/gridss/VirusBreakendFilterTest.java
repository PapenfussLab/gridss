package gridss;

import au.edu.wehi.idsv.Hg19Tests;
import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class VirusBreakendFilterTest extends IntermediateFilesTest {
    @Test
    public void should_filter_no_host_alignment() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("no_host_alignment_host")));
    }
    @Test
    public void should_filter_repeat_50_percent_overlap() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.MINIMUM_REPEAT_OVERLAP = 0.49;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("repeat_overlap_host")));
    }
    @Test
    public void should_filter_short_host_hits() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("short_host_alignment_host")));
    }
    @Test
    public void should_filter_breakpoints() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("breakpoint_host")));
    }
    @Test
    public void should_retain_small_contained_repeat() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertTrue(vcf.stream().anyMatch(vc -> vc.getID().equals("small_contained_repeat_host")));
    }
    @Test
    public void should_retain_simple_repeat_no_overlap() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertTrue(vcf.stream().anyMatch(vc -> vc.getID().equals("repeat_no_overlap_host")));
    }
    @Test
    public void should_retain_complex_repeat_overlap() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertTrue(vcf.stream().anyMatch(vc -> vc.getID().equals("line_overlap_host")));
    }
    @Test
    public void insert_sequence_should_match() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/test.vcf");
        cmd.OUTPUT = output;
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();

        List<VariantContext> vcf = getRawVcf(output);
        VariantContext vc = vcf.stream().filter(x -> x.getID().equals("virus_forward_host_plus_virus")).findFirst().get();
        //kraken:taxid|333760|NC_001526	1037	virus_forward_host_plus	T	NACGT.	100	PASS	BEALN=polyA:10|+|1S3M|60
        Assert.assertEquals("NA[polyA:10[", vc.getAlternateAllele(0).getDisplayString());
        vc = vcf.stream().filter(x -> x.getID().equals("virus_forward_host_plus_host")).findFirst().get();
        Assert.assertEquals("]kraken:taxid|333760|NC_001526:1037]AN", vc.getAlternateAllele(0).getDisplayString());
        Assert.assertEquals("polyA", vc.getContig());
        Assert.assertEquals(10, vc.getStart());

        //kraken:taxid|333760|NC_001526	1037	virus_forward_host_minus	T	NACGT.	100	PASS	BEALN=polyA:10|-|3M1S|60
        vc = vcf.stream().filter(x -> x.getID().equals("virus_forward_host_minus_virus")).findFirst().get();
        Assert.assertEquals("NA]polyA:12]", vc.getAlternateAllele(0).getDisplayString());
        vc = vcf.stream().filter(x -> x.getID().equals("virus_forward_host_minus_host")).findFirst().get();
        Assert.assertEquals("NT]kraken:taxid|333760|NC_001526:1037]", vc.getAlternateAllele(0).getDisplayString());

        //kraken:taxid|333760|NC_001526	1037	virus_back_host_plus	T	.AACTGG	100	PASS	BEALN=polyA:10|+|3M2S|60
        vc = vcf.stream().filter(x -> x.getID().equals("virus_back_host_plus_virus")).findFirst().get();
        Assert.assertEquals("]polyA:12]TGG", vc.getAlternateAllele(0).getDisplayString());
        vc = vcf.stream().filter(x -> x.getID().equals("virus_back_host_plus_host")).findFirst().get();
        Assert.assertEquals("NTG[kraken:taxid|333760|NC_001526:1037[", vc.getAlternateAllele(0).getDisplayString());

        // kraken:taxid|333760|NC_001526	1037	virus_back_host_minus	T	.AACTGG	100	PASS	BEALN=polyA:10|-|2S3M|60
        vc = vcf.stream().filter(x -> x.getID().equals("virus_back_host_minus_virus")).findFirst().get();
        Assert.assertEquals("[polyA:10[TGG", vc.getAlternateAllele(0).getDisplayString());
        vc = vcf.stream().filter(x -> x.getID().equals("virus_back_host_minus_host")).findFirst().get();
        Assert.assertEquals("[kraken:taxid|333760|NC_001526:1037[CAN", vc.getAlternateAllele(0).getDisplayString());
    }
    @Test
    public void should_match_host_taxid() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.TAXONOMY_IDS = ImmutableList.of(9606);
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertTrue(vcf.stream().anyMatch(vc -> vc.getID().equals("correct_host_taxid_host")));
        Assert.assertFalse(vcf.stream().anyMatch(vc -> vc.getID().equals("wrong_host_taxid_host")));
    }
    @Test
    public void should_apply_LOW_MAPQ_filter() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/candidates.vcf");
        cmd.OUTPUT = output;
        cmd.TAXONOMY_IDS = ImmutableList.of(9606);
        cmd.setReference(SMALL_FA_FILE);
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
        Assert.assertFalse(vcf.stream().filter(vc -> vc.getID().equals("small_contained_repeat_host")).findFirst().get().getFilters().contains("LOW_MAPQ"));
        Assert.assertFalse(vcf.stream().filter(vc -> vc.getID().equals("small_contained_repeat_virus")).findFirst().get().getFilters().contains("LOW_MAPQ"));
        Assert.assertTrue(vcf.stream().filter(vc -> vc.getID().equals("correct_host_taxid_host")).findFirst().get().getFilters().contains("LOW_MAPQ"));
        Assert.assertTrue(vcf.stream().filter(vc -> vc.getID().equals("correct_host_taxid_virus")).findFirst().get().getFilters().contains("LOW_MAPQ"));
    }
    @Test
    @Category(Hg19Tests.class)
    public void regression_should_not_crash_when_viral_contig_contains_period() {
        VirusBreakendFilter cmd = new VirusBreakendFilter();
        File output = new File(testFolder.getRoot(), "out.vcf");
        cmd.INPUT = new File("src/test/resources/virusbreakend/VirusBreakendFilter_regression.vcf");
        cmd.OUTPUT = output;
        cmd.TAXONOMY_IDS = ImmutableList.of(9606);
        cmd.setReference(Hg19Tests.findBroadHg19Reference());
        cmd.doWork();
        List<VariantContext> vcf = getRawVcf(output);
    }
}