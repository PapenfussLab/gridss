package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Test;

public class SocratesVariantContextTest extends TestHelper {
	@Test
	public void evidence_should_start_null() {
		SocratesVariantContext vc = SocratesVariantContext.create(null, minimalVariant().make());
		assertNull(vc.getReferenceSpanningPairs());
		assertNull(vc.getReferenceReadDepth());
		assertNull(vc.getOEACount());
		assertNull(vc.getSoftClipCount());
		assertNull(vc.getDiscordantPairCount());
	}
	@Test
	public void addEvidenceAttributes_should_set_evidence() {
		SocratesVariantContext vc = SocratesVariantContext.create(null, minimalVariant().make());
		vc = vc.addEvidenceAttributes(vc, 1, 2, 3, 4, 5);
		assertEquals(1, (int)vc.getReferenceSpanningPairs());
		assertEquals(2, (int)vc.getReferenceReadDepth());
		assertEquals(3, (int)vc.getOEACount());
		assertEquals(4, (int)vc.getSoftClipCount());
		assertEquals(5, (int)vc.getDiscordantPairCount());
	}
	@Test
	public void getReferenceIndex_should_lookup_dictionary_of_chr() {
		VariantContext vc = minimalVariant()
			.chr("polyA")
			.make();
		assertEquals(0, SocratesVariantContext.create(getContext(), vc).getReferenceIndex());
	}
	@Test
	public void create_should_default_to_SocratesVariantContext() {
		assertTrue(SocratesVariantContext.create(getContext(), minimalVariant().make()) instanceof SocratesVariantContext);
	}
	@Test
	public void create_should_make_DirectedBreakpointAssembly() {
		DirectedBreakpointAssembly vc = DirectedBreakpointAssembly.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("C"), B("AC"), 1, 0);
		assertTrue(SocratesVariantContext.create(getContext(), vc) instanceof DirectedBreakpointAssembly);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint() {
		VariantContext vc = minimalVariant().alleles("A", "A[polyA:1[").make();
		assertTrue(SocratesVariantContext.create(getContext(), vc) instanceof VariantContextDirectedBreakpoint);
	}
}
