package au.edu.wehi.idsv.configuration;

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
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 100, 200), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 2), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 3, 3), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, FWD, 2, 2), "")).contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void breakpointFilters_should_filter_insertion() {
		VariantCallingConfiguration p = getConfig().getVariantCalling();
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 100, 200), "A")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 2), "A")).contains(VcfFilter.REFERENCE_ALLELE));
	}
}
