package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.StructuralVariationCallBuilder;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.vcf.VcfAttributes;



public class StructuralVariationCallBuilderTest extends TestHelper {
	@Test
	public void should_sum_evidence_metrics() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), new BreakendSummary(0, FWD, 1, 1, null));
		builder.evidence(SCE(FWD, Read(0, 1, "1M5S")));
		builder.evidence(SCE(FWD, Read(0, 1, "1M6S")));
		VariantContextDirectedBreakpoint bp = builder.make();
		assertEquals(2, bp.getAttributeAsInt(VcfAttributes.SOFT_CLIP_READ_COUNT.attribute(), 0));
		assertEquals(11, bp.getAttributeAsInt(VcfAttributes.SOFT_CLIP_TOTAL_LENGTH.attribute(), 0));
	}
	@Test(expected=IllegalArgumentException.class)
	public void evidence_must_support_call() {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(getContext(), new BreakendSummary(0, FWD, 1, 1, null));
		builder.evidence(SCE(FWD, Read(0, 2, "1M5S")));
		builder.make();
	}
}
