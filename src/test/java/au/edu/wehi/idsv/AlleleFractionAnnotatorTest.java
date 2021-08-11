package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Assert;
import org.junit.Test;

public class AlleleFractionAnnotatorTest extends TestHelper {
    @Test
    public void should_include_rp_for_large_events() {
        ProcessingContext pc = getContext();
        IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext(), new StructuralVariationCallBuilder(pc, new CalledBreakpointPositionLookup(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
            breakpoint(new BreakpointSummary(2, FWD, 78, 6, BWD, 79), "");
            phredScore(50);
        }}.make()).make());
        builder.genotypeBuilder.get(0).attribute("REF", 2);
        builder.genotypeBuilder.get(0).attribute("REFPAIR", 2);
        builder.genotypeBuilder.get(0).attribute("VF", 4);
        builder.genotypeBuilder.get(1).attribute("REF", 3);
        builder.genotypeBuilder.get(1).attribute("REFPAIR", 3);
        builder.genotypeBuilder.get(1).attribute("VF", 2);
        VariantContextDirectedBreakpoint e = (VariantContextDirectedBreakpoint)builder.make();
        MockSAMEvidenceSource cat1 = SES(0, 2000);
        cat1.category = 1;
        AlleleFractionAnnotator aaf = new AlleleFractionAnnotator(getContext(), ImmutableList.of(SES(0, 500), cat1));
        VariantContext var = aaf.annotate(e);
        Assert.assertEquals(0.5, (double)var.getGenotype(0).getAnyAttribute("AF"), 0);
        Assert.assertEquals(0.25, (double)var.getGenotype(1).getAnyAttribute("AF"), 0);
    }
    @Test
    public void should_exclude_rp_for_small_events() {
        ProcessingContext pc = getContext();
        IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext(), new StructuralVariationCallBuilder(pc, new CalledBreakpointPositionLookup(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
            breakpoint(new BreakpointSummary(0, FWD, 100, 0, BWD, 500+200-1), "");
            phredScore(50);
        }}.make()).make());
        builder.genotypeBuilder.get(0).attribute("REF", 4);
        builder.genotypeBuilder.get(0).attribute("REFPAIR", 5);
        builder.genotypeBuilder.get(0).attribute("VF", 1);
        builder.genotypeBuilder.get(1).attribute("REF", 2);
        builder.genotypeBuilder.get(1).attribute("REFPAIR", 3);
        builder.genotypeBuilder.get(1).attribute("VF", 6);
        VariantContextDirectedBreakpoint e = (VariantContextDirectedBreakpoint)builder.make();
        MockSAMEvidenceSource cat1 = SES(200, 400);
        cat1.category = 1;
        AlleleFractionAnnotator aaf = new AlleleFractionAnnotator(getContext(), ImmutableList.of(SES(400, 500), cat1));
        VariantContext var = aaf.annotate(e);
        Assert.assertEquals((double)var.getGenotype(0).getAnyAttribute("AF"), 0.2,0);
        Assert.assertEquals((double)var.getGenotype(1).getAnyAttribute("AF"), 0.75, 0);
    }
    @Test
    public void should_use_BVF_for_breakends() {
        ProcessingContext pc = getContext();
        IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext(), new StructuralVariationCallBuilder(pc, new CalledBreakpointPositionLookup(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
            breakend(new BreakendSummary(0, FWD, 100), "ACGT");
            phredScore(50);
        }}.make()).make());
        builder.genotypeBuilder.get(0).attribute("REF", 4);
        builder.genotypeBuilder.get(0).attribute("REFPAIR", 5);
        builder.genotypeBuilder.get(0).attribute("BVF", 1);
        builder.genotypeBuilder.get(1).attribute("REF", 2);
        builder.genotypeBuilder.get(1).attribute("REFPAIR", 2);
        builder.genotypeBuilder.get(1).attribute("BVF", 6);
        VariantContextDirectedEvidence e = (VariantContextDirectedEvidence) builder.make();
        MockSAMEvidenceSource cat1 = SES(0, 2000);
        cat1.category = 1;
        AlleleFractionAnnotator aaf = new AlleleFractionAnnotator(getContext(), ImmutableList.of(SES(0, 500), cat1));
        VariantContext var = aaf.annotate(e);
        Assert.assertEquals((double)var.getGenotype(0).getAnyAttribute("AF"), 0.1,0);
        Assert.assertEquals((double)var.getGenotype(1).getAnyAttribute("AF"), 0.6, 0);
    }
}