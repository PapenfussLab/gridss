package au.edu.wehi.idsv.debruijn;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;

public class ContigCategorySupportHelperTest extends TestHelper {
	@Test
	public void should_convert_to_match_mismatch_cigar() {
		List<BitSet> list = new ArrayList<>();
		list.add(new BitSet(4));
		list.get(0).set(0);
		list.get(0).set(2);
		list.add(new BitSet(4));
		list.get(1).set(1);
		list.get(1).set(2);
		list.add(new BitSet(4));
		String cigar = ContigCategorySupportHelper.asSupportCigars(4, list);
		Assert.assertEquals("1=1X1=1X,1X2=1X,4X", cigar);
	}
	@Test
	public void should_decode_cigar() {
		List<BitSet> list = ContigCategorySupportHelper.asPerCategoryReadCoverage("1=1X1=1X,1X2=1X,4X");
		Assert.assertEquals(3, list.size());
		Assert.assertTrue(list.get(0).get(0));
		Assert.assertTrue(list.get(0).get(2));
		Assert.assertTrue(list.get(1).get(1));
		Assert.assertTrue(list.get(1).get(2));
		Assert.assertEquals(2, list.get(0).cardinality());
		Assert.assertEquals(2, list.get(1).cardinality());
		Assert.assertEquals(0, list.get(2).cardinality());
	}
	@Test
	public void should_require_coverage_on_both_sides_of_breakend_position() {
		Assert.assertEquals(ImmutableList.of(false, false), ContigCategorySupportHelper.supportsBreakendBefore(0, "4X,4X"));
		Assert.assertEquals(ImmutableList.of(false, true), ContigCategorySupportHelper.supportsBreakendBefore(1, "4X,4="));
		Assert.assertEquals(ImmutableList.of(false, false), ContigCategorySupportHelper.supportsBreakendBefore(1, "4X,1X1=2X"));
		Assert.assertEquals(ImmutableList.of(false, true), ContigCategorySupportHelper.supportsBreakendBefore(1, "4X,2=2X"));
	}
	@Test
	public void supportsBreakendBefore_should_not_support_out_of_bounds_positions() {
		Assert.assertEquals(ImmutableList.of(false), ContigCategorySupportHelper.supportsBreakendBefore(-1, "4="));
		Assert.assertEquals(ImmutableList.of(false), ContigCategorySupportHelper.supportsBreakendBefore(0, "4="));
		Assert.assertEquals(ImmutableList.of(true), ContigCategorySupportHelper.supportsBreakendBefore(1, "4="));
		Assert.assertEquals(ImmutableList.of(true), ContigCategorySupportHelper.supportsBreakendBefore(2, "4="));
		Assert.assertEquals(ImmutableList.of(true), ContigCategorySupportHelper.supportsBreakendBefore(3, "4="));
		Assert.assertEquals(ImmutableList.of(false), ContigCategorySupportHelper.supportsBreakendBefore(4, "4="));
		Assert.assertEquals(ImmutableList.of(false), ContigCategorySupportHelper.supportsBreakendBefore(5, "4="));
	}
}
