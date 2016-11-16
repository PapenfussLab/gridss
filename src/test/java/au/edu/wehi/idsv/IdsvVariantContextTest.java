package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class IdsvVariantContextTest extends TestHelper {
	@Test
	public void getReferenceIndex_should_lookup_dictionary_of_chr() {
		VariantContext vc = minimalVariant()
			.chr("polyA")
			.make();
		assertEquals(0, IdsvVariantContext.create(getContext(), AES(), vc).getReferenceIndex());
	}
	@Test
	public void create_should_default_to_IdsvVariantContext() {
		assertTrue(IdsvVariantContext.create(getContext(), AES(), minimalVariant().make()) instanceof IdsvVariantContext);
	}
	@Test
	public void create_should_make_VariantContextDirectedEvidence_from_breakend() {
		VariantContext vc = minimalVariant().alleles("A", "A.").make();
		assertTrue(IdsvVariantContext.create(getContext(), AES(), vc) instanceof VariantContextDirectedEvidence);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint_from_breakpoint() {
		VariantContext vc = minimalVariant().alleles("A", "A[polyA:1[").make();
		assertTrue(IdsvVariantContext.create(getContext(), AES(), vc) instanceof VariantContextDirectedBreakpoint);
	}
	public class TestIdsvVariantContext extends IdsvVariantContext {
		public TestIdsvVariantContext(VariantContext context) {
			super(getContext(), AES(), context);
		}
		public TestIdsvVariantContext(GenomicProcessingContext processContext, VariantContext context) {
			super(processContext, AES(), context);
		}
		public List<Integer> getAttributeAsIntList(String attrName) {
			return AttributeConverter.asIntList(getAttribute(attrName));
		}
		public List<String> getAttributeAsStringList(String attrName) {
			return AttributeConverter.asStringList(getAttribute(attrName));
		}
		public int getAttributeAsIntListOffset(String attrName, int i, int j) {
			return AttributeConverter.asIntListOffset(getAttribute(attrName), i, j);
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
	@Test
	public void getAttributeAsStringList() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().attribute("slist", L("1", "2")).make());
		assertEquals(2, vc.getAttributeAsStringList("slist").size());
	}
	@Test
	public void getAttributeAsIntListOffset_should_default_if_not_enough_elements() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().attribute("intlist", L("1", "2")).make());
		assertEquals(7, vc.getAttributeAsIntListOffset("intlist", 3, 7));
	}
	@Test
	public void getAttributeAsIntListOffset_should_default_if_not_attribute() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().make());
		assertEquals(7, vc.getAttributeAsIntListOffset("intlist", 3, 7));
	}
	@Test
	public void getAttributeAsIntListOffset_should_default_should_get_ith_element() {
		TestIdsvVariantContext vc = new TestIdsvVariantContext(minimalVariant().attribute("intlist", L("1", "2")).make());
		assertEquals(1, vc.getAttributeAsIntListOffset("intlist", 0, 7));
	}
}
