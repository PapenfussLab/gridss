package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

import org.junit.Test;

public class IdsvVariantContextTest extends TestHelper {
	@Test
	public void getReferenceIndex_should_lookup_dictionary_of_chr() {
		VariantContext vc = minimalVariant()
			.chr("polyA")
			.make();
		assertEquals(0, IdsvVariantContext.create(getContext(), AES(), vc).getReferenceIndex());
	}
	@Test
	public void create_should_default_to_SocratesVariantContext() {
		assertTrue(IdsvVariantContext.create(getContext(), AES(), minimalVariant().make()) instanceof IdsvVariantContext);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint() {
		VariantContextDirectedEvidence vc = AE();
		assertTrue(IdsvVariantContext.create(getContext(), AES(), vc) instanceof VariantContextDirectedEvidence);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint_from_breakpoint() {
		VariantContext vc = minimalVariant().alleles("A", "A[polyA:1[").make();
		assertTrue(IdsvVariantContext.create(getContext(), AES(), vc) instanceof VariantContextDirectedEvidence);
	}
	public class TestIdsvVariantContext extends IdsvVariantContext {
		public TestIdsvVariantContext(VariantContext context) {
			super(getContext(), AES(), context);
		}
		public TestIdsvVariantContext(ProcessingContext processContext, VariantContext context) {
			super(processContext, AES(), context);
		}
		@Override
		public List<Integer> getAttributeAsIntList(String attrName) {
			return super.getAttributeAsIntList(attrName);
		}
	}
	@Test
	public void getAttributeAsIntList_should_allow_array() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().attribute("intlist", new int[] { 1, 2}).make());
		assertEquals(2, vc.getAttributeAsIntList("intlist").size());
	}
	@Test
	public void getAttributeAsIntList_should_allow_int_list() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().attribute("intlist", L(1, 2)).make());
		assertEquals(2, vc.getAttributeAsIntList("intlist").size());
	}
	@Test
	public void getAttributeAsIntList_should_allow_single_int() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().attribute("intlist", 1).make());
		assertEquals(1, vc.getAttributeAsIntList("intlist").size());
	}
	/**
	 * Picard parses into ArrayList of String
	 */
	@Test
	public void getAttributeAsIntList_should_allow_string_list() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().attribute("intlist", L("1", "2")).make());
		assertEquals(2, vc.getAttributeAsIntList("intlist").size());
	}
}
