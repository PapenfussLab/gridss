package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;


public class VariantContextDirectedBreakpointTest extends TestHelper {
	public VariantContextDirectedBreakpoint BP(String id, BreakpointSummary bp) {
		IdsvVariantContextBuilder builder = minimalBreakend()
				.breakpoint(bp, "");
		builder.id(id);
		return (VariantContextDirectedBreakpoint)builder.make();
	}
	@Test
	public void ByRemoteBreakendLocationStart_should_sort_by_remote_breakend() {
		List<VariantContextDirectedBreakpoint> list = Lists.newArrayList(
				BP("1", new BreakpointSummary(0, FWD, 1, 1, 2, FWD, 10, 10)),
				BP("2", new BreakpointSummary(1, FWD, 1, 1, 1, FWD, 10, 10)),
				BP("3", new BreakpointSummary(1, FWD, 5, 5, 1, FWD, 1, 1)),
				BP("4", new BreakpointSummary(1, FWD, 10, 10, 0, FWD, 1, 1)),
				BP("5", new BreakpointSummary(2, FWD, 1, 1, 1, FWD, 2, 2))
				);
		Collections.sort(list, VariantContextDirectedBreakpoint.ByRemoteBreakendLocationStart);
		List<VariantContextDirectedBreakpoint> result = list;
		
		assertEquals(5, result.size());
		assertEquals("4", result.get(0).getID());
		assertEquals("3", result.get(1).getID());
		assertEquals("5", result.get(2).getID());
		assertEquals("2", result.get(3).getID());
		assertEquals("1", result.get(4).getID());
	}
	@Test
	public void ByRemoteBreakendLocationStartRaw_should_sort_by_remote_breakend() {
		List<VariantContext> list = Lists.newArrayList(
				(VariantContext)BP("1", new BreakpointSummary(0, FWD, 1, 1, 2, FWD, 10, 10)),
				(VariantContext)BP("2", new BreakpointSummary(1, FWD, 1, 1, 1, FWD, 10, 10)),
				(VariantContext)BP("3", new BreakpointSummary(1, FWD, 5, 5, 1, FWD, 1, 1)),
				(VariantContext)BP("4", new BreakpointSummary(1, FWD, 10, 10, 0, FWD, 1, 1)),
				(VariantContext)BP("5", new BreakpointSummary(2, FWD, 1, 1, 1, FWD, 2, 2))
				);
		Collections.sort(list, VariantContextDirectedBreakpoint.ByRemoteBreakendLocationStartRaw(getContext()));
		List<VariantContext> result = list;
		
		assertEquals(5, result.size());
		assertEquals("4", result.get(0).getID());
		assertEquals("3", result.get(1).getID());
		assertEquals("5", result.get(2).getID());
		assertEquals("2", result.get(3).getID());
		assertEquals("1", result.get(4).getID());
	}
}
