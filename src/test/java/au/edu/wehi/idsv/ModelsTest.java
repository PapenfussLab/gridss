package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;

import org.junit.Test;

import com.google.common.collect.Lists;

public class ModelsTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void calculateBreakend_should_throw_upon_null_evidence() {
		Models.calculateBreakend(getContext().getLinear(), null);
	}
	@Test(expected=IllegalArgumentException.class)
	public void calculateBreakend_should_throw_upon_no_evidence() {
		Models.calculateBreakend(getContext().getLinear(), new ArrayList<DirectedEvidence>());
	}
	@Test
	public void calculateBreakend_should_return_breakend() {
		BreakendSummary bs = Models.calculateBreakend(getContext().getLinear(), Lists.newArrayList(
			new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 1, 2, 1, BWD, 2, 2))
				));
		assertEquals(new BreakendSummary(0, FWD, 1, 2), bs);
		assertEquals(BreakendSummary.class, bs.getClass());
	}
	@Test
	public void calculateBreakend_should_return_interval_for_single_evidence() {
		BreakendSummary bs = Models.calculateBreakend(getContext().getLinear(), Lists.newArrayList(
			new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 1, 2, 1, BWD, 2, 2))
				));
		assertEquals(new BreakendSummary(0, FWD, 1, 2), bs);
	}
	@Test
	public void calculateBreakend_should_return_overlap_for_multiple_evidence() {
		BreakendSummary bs = Models.calculateBreakend(getContext().getLinear(), Lists.newArrayList(
			new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 1, 4, 1, BWD, 2, 2)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 3, 5))
				));
		assertEquals(new BreakendSummary(0, FWD, 3, 4), bs);
	}
	@Test
	public void calculateBreakend_should_ignore_overlaps_resulting_in_no_breakend() {
		BreakendSummary bs = Models.calculateBreakend(getContext().getLinear(), Lists.newArrayList(
			new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 1, 4, 1, BWD, 2, 2)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 3, 5)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 5, 5))
				));
		assertEquals(new BreakendSummary(0, FWD, 3, 4), bs);
	}
	@Test
	public void calculateBreakend_should_reduce_to_consistent_set() {
		BreakendSummary bs = Models.calculateBreakend(getContext().getLinear(), Lists.newArrayList(
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 10, 15)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 12, 18))
				));
		assertEquals(new BreakendSummary(0, FWD, 12, 15), bs);
	}
	@Test
	public void calculateBreakend_calculate_interval_with_greatest_support() {
		//          1         2
		// 12345678901234567890123456789
		//          ****** **@@
		//               ****@@@
		//         *******   
		//               ^
		//            most support
		BreakendSummary bs = Models.calculateBreakend(getContext().getLinear(), Lists.newArrayList(
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 10, 15)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 15, 18)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 17, 18)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 9, 15)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 19, 20)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 19, 21))
				));
		assertEquals(new BreakendSummary(0, FWD, 15, 15), bs);
	}
	@Test
	public void calculateBreakend_should_infer_direction_from_support() {
		BreakendSummary bs = Models.calculateBreakend(getContext().getLinear(), Lists.newArrayList(
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 10, 15)),
			new MockDirectedEvidence(new BreakendSummary(0, FWD, 10, 16)),
			new MockDirectedEvidence(new BreakendSummary(0, BWD, 12, 18))
				));
		assertEquals(new BreakendSummary(0, FWD, 10, 15), bs);
	}
}
