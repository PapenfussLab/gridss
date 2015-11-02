package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import htsjdk.variant.vcf.VCFConstants;


public class AttributeConverterTest {
	@Test
	public void asIntList_should_allow_array() {
		assertEquals(2, AttributeConverter.asIntList(new int[] { 1, 2}, null).size());
	}
	@Test
	public void asIntList_should_allow_int_list() {
		assertEquals(2, AttributeConverter.asIntList(ImmutableList.<Integer>of(1, 2), null).size());
	}
	@Test
	public void asIntList_should_allow_single_int() {
		assertEquals(1, AttributeConverter.asIntList(7, null).size());
		assertEquals(7, (int)AttributeConverter.asIntList(7, null).get(0));
	}
	/**
	 * Picard parses into ArrayList of String
	 */
	@Test
	public void asIntList_should_allow_string_list() {
		assertEquals(2, AttributeConverter.asIntList(ImmutableList.<String>of("1", "2"), null).size());
	}
	@Test
	public void asStringList() {
		assertEquals(2, AttributeConverter.asStringList(ImmutableList.<String>of("1", "2"), null).size());
	}
	@Test
	public void asIntListOffset_should_get_offset() {
		assertEquals(2, AttributeConverter.asIntListOffset(ImmutableList.<Integer>of(1, 2, 7), 1, 0));
	}
	@Test
	public void asIntListOffset_should_default_if_not_enough_elements() {
		assertEquals(7, AttributeConverter.asIntListOffset(ImmutableList.<Integer>of(1, 2), 3, 7));
	}
	@Test
	public void asIntListOffset_should_if_null() {
		assertEquals(7, AttributeConverter.asIntListOffset(null, 3, 7));
	}
	@Test
	public void vcf_null_should_be_treated_as_null() {
		assertEquals("default", AttributeConverter.asStringList(VCFConstants.MISSING_VALUE_v4, ImmutableList.of("default")).get(0));
	}
	@Test
	public void asInt_should_return_first_value_of_list() {
		assertEquals(Integer.valueOf(7), AttributeConverter.asInt(ImmutableList.<Object>of("7", 2), 0));
	}
}
