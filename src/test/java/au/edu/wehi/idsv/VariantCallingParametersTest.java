package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class VariantCallingParametersTest extends TestHelper {
	@Test
	public void breakpointFilters_should_filter_small_indels() {
		VariantCallingParameters p = new VariantCallingParameters();
		p.minIndelSize = 3;
		assertEquals(0, p.breakpointFilters(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 1, 1)).size()); // not an indel
		assertEquals(0, p.breakpointFilters(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 2)).size()); // indel size of 0 - ie: no variant
		assertEquals(1, p.breakpointFilters(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 3, 3)).size()); // 1
		assertEquals(1, p.breakpointFilters(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 4, 4)).size()); // 2
		assertEquals(0, p.breakpointFilters(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 5, 5)).size()); // 3 min size is not filtered, only less than min
	}
	@Test
	public void breakpointFilters_should_filter_if_possibly_small_indel() {
		VariantCallingParameters p = new VariantCallingParameters();
		p.minIndelSize = 10;
		assertEquals(1, p.breakpointFilters(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 1, 102)).size()); // possibly 1bp del of base 101
		assertEquals(1, p.breakpointFilters(new BreakpointSummary(0, FWD, 1, 100, 0, BWD, 110, 200)).size()); // possibly 9bp
		assertEquals(1, p.breakpointFilters(new BreakpointSummary(0, FWD, 100, 400, 0, BWD, 200, 200)).size());
		assertEquals(0, p.breakpointFilters(new BreakpointSummary(0, FWD, 100, 200, 0, BWD, 1, 101)).size()); // possibly 0bp (ie: no variant)
	}
}
