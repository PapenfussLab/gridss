package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class IdsvVariantContextTest extends TestHelper {
	@Test
	public void getReferenceIndex_should_lookup_dictionary_of_chr() {
		VariantContext vc = minimalVariant()
			.chr("polyA")
			.make();
		assertEquals(0, IdsvVariantContext.create(getContext().getDictionary(), AES(), vc).getReferenceIndex());
	}
	@Test
	public void create_should_default_to_IdsvVariantContext() {
		assertTrue(IdsvVariantContext.create(getContext().getDictionary(), AES(), minimalVariant().make()) instanceof IdsvVariantContext);
	}
	@Test
	public void create_should_make_VariantContextDirectedEvidence_from_breakend() {
		VariantContext vc = minimalVariant().alleles("A", "A.").make();
		assertTrue(IdsvVariantContext.create(getContext().getDictionary(), AES(), vc) instanceof VariantContextDirectedEvidence);
	}
	@Test
	public void create_should_make_VariantContextDirectedBreakpoint_from_breakpoint() {
		VariantContext vc = minimalVariant().alleles("A", "A[polyA:1[").make();
		assertTrue(IdsvVariantContext.create(getContext().getDictionary(), AES(), vc) instanceof VariantContextDirectedBreakpoint);
	}
	public class TestIdsvVariantContext extends IdsvVariantContext {
		private static final long serialVersionUID = 1L;
		public TestIdsvVariantContext(VariantContext context) {
			super(getContext().getDictionary(), AES(), context);
		}
		public TestIdsvVariantContext(GenomicProcessingContext processContext, VariantContext context) {
			super(processContext.getDictionary(), AES(), context);
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
	@Test
	public void ByLocation_should_have_stable_sort_order() {
		IdsvVariantContext v1 = (IdsvVariantContext)new IdsvVariantContextBuilder(getContext()).chr("polyA").start(1).stop(1).alleles("A", "C").make();
		IdsvVariantContext v2 = (IdsvVariantContext)new IdsvVariantContextBuilder(getContext()).chr("polyA").start(1).stop(1).alleles("A", "T").make();
		for (Ordering<IdsvVariantContext> order : ImmutableList.of(IdsvVariantContext.ByLocationStart, IdsvVariantContext.ByLocationEnd)) {
			List<IdsvVariantContext> list1 = new ArrayList<>();
			list1.add(v1);
			list1.add(v2);
			List<IdsvVariantContext> list2 = new ArrayList<>();
			list2.add(v2);
			list2.add(v1);
			list1.sort(order);
			list2.sort(order);
			assertEquals(list1, list2);
		}
	}
}
