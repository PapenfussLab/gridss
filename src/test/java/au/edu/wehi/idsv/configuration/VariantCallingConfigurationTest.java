package au.edu.wehi.idsv.configuration;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.vcf.VcfFilter;


public class VariantCallingConfigurationTest extends TestHelper {
	public VariantContextDirectedBreakpoint V(BreakpointSummary bs, String untemplated) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.breakpoint(bs, untemplated);
		return (VariantContextDirectedBreakpoint)builder.make();
	}
	@Test
	public void breakpointFilters_should_filter_reference_allele() {
		VariantCallingConfiguration p = getConfig().getVariantCalling();
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 150, 100, 200, 0, BWD, 150, 100, 200), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 1, 0, BWD, 2, 2, 2), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 1, 0, BWD, 3, 3, 3), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 1, 0, FWD, 2, 2, 2), "")).contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void breakpointFilters_should_filter_insertion() {
		VariantCallingConfiguration p = getConfig().getVariantCalling();
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 150, 100, 200, 0, BWD, 150, 100, 200), "A")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 1, 0, BWD, 2, 2, 2), "A")).contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void full_margin_should_apply_to_interchromosomal_events() {
		VariantCallingConfiguration config = new VariantCallingConfiguration(getDefaultConfig());
		config.breakendMargin = 100;
		assertEquals(new BreakpointSummary(0, FWD, 1000, 900, 1100, 1, FWD, 1000, 900, 1100), config.withMargin(new BreakpointSummary(0, FWD, 1000, 1, FWD, 1000)));
	}
	@Test
	public void margin_should_linearly_reduce_from_margin_at_twice_margin_size_to_zero_at_1bp() {
		VariantCallingConfiguration config = new VariantCallingConfiguration(getDefaultConfig());
		config.breakendMargin = 10;
		int[] size = new int[] { 20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0 };
		int[] expectedMargin = new int[] { 10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0, };
		for (int i = 0; i < size.length; i++) {
			assertEquals(new BreakpointSummary(0, FWD, 100, 100-expectedMargin[i], 100+expectedMargin[i], 0, FWD, 100 + size[i], 100-expectedMargin[i] + size[i], 100+expectedMargin[i] + size[i]),
					config.withMargin(new BreakpointSummary(0, FWD, 100, 100, 100, 0, FWD, 100 + size[i], 100 + size[i], 100 + size[i])));
		}
	}
}
