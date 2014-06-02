package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.variant.variantcontext.VariantContext;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ReadEvidenceAssemblerUtil;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;

public class IdsvVariantContextTest extends TestHelper {
	@Test
	public void getReferenceIndex_should_lookup_dictionary_of_chr() {
		VariantContext vc = minimalVariant()
			.chr("polyA")
			.make();
		assertEquals(0, IdsvVariantContext.create(getContext(), vc).getReferenceIndex());
	}
	@Test
	public void create_should_default_to_SocratesVariantContext() {
		assertTrue(IdsvVariantContext.create(getContext(), minimalVariant().make()) instanceof IdsvVariantContext);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint() {
		VariantContextDirectedBreakpoint vc = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("C"), B("AC"), 1, 1, 0);
		assertTrue(IdsvVariantContext.create(getContext(), vc) instanceof VariantContextDirectedBreakpoint);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint_from_breakpoint() {
		VariantContext vc = minimalVariant().alleles("A", "A[polyA:1[").make();
		assertTrue(IdsvVariantContext.create(getContext(), vc) instanceof VariantContextDirectedBreakpoint);
	}
}
