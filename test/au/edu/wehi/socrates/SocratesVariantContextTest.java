package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Test;

public class SocratesVariantContextTest extends TestHelper {
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
	public void create_should_make_VariantContextDirectedBreakpoint() {
		VariantContextDirectedBreakpoint vc = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("C"), B("AC"), 1, 1, 0);
		assertTrue(SocratesVariantContext.create(getContext(), vc) instanceof VariantContextDirectedBreakpoint);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint_from_breakpoint() {
		VariantContext vc = minimalVariant().alleles("A", "A[polyA:1[").make();
		assertTrue(SocratesVariantContext.create(getContext(), vc) instanceof VariantContextDirectedBreakpoint);
	}
}
