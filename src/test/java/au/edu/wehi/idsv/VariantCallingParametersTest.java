package au.edu.wehi.idsv;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfFilter;


public class VariantCallingParametersTest extends TestHelper {
	public VariantContextDirectedBreakpoint V(BreakpointSummary bs, String untemplated) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder.breakpoint(bs, untemplated);
		return (VariantContextDirectedBreakpoint)builder.make();
	}
	@Test
	public void breakpointFilters_should_filter_small_indels() {
		VariantCallingParameters p = new VariantCallingParameters();
		p.minIndelSize = 3;
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 1, 1), "")).contains(VcfFilter.SMALL_INDEL)); // not an indel
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 2), "")).contains(VcfFilter.SMALL_INDEL)); // indel size of 0 - ie: no variant
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 3, 3), "")).contains(VcfFilter.SMALL_INDEL)); // 1
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 4, 4), "")).contains(VcfFilter.SMALL_INDEL)); // 2
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 5, 5), "")).contains(VcfFilter.SMALL_INDEL)); // 3 min size is not filtered, only less than min
	}
	@Test
	public void breakpointFilters_should_filter_if_possibly_small_indel() {
		VariantCallingParameters p = new VariantCallingParameters();
		p.minIndelSize = 10;
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 1, 102), "")).contains(VcfFilter.SMALL_INDEL)); // possibly 1bp del of base 101
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 100, 0, BWD, 110, 200), "")).contains(VcfFilter.SMALL_INDEL)); // possibly 9bp
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 100, 400, 0, BWD, 200, 200), "")).contains(VcfFilter.SMALL_INDEL));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 1, 101), "")).contains(VcfFilter.SMALL_INDEL)); // possibly 0bp (ie: no variant)
	}
	@Test
	public void breakpointFilters_should_filter_reference_allele() {
		VariantCallingParameters p = new VariantCallingParameters();
		p.minIndelSize = 10;
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 100, 200), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertTrue(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 2), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 3, 3), "")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, FWD, 2, 2), "")).contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void breakpointFilters_should_filter_insertion() {
		VariantCallingParameters p = new VariantCallingParameters();
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 100, 200), "A")).contains(VcfFilter.REFERENCE_ALLELE));
		assertFalse(p.calculateBreakpointFilters(V(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 2), "A")).contains(VcfFilter.REFERENCE_ALLELE));
	}
}
