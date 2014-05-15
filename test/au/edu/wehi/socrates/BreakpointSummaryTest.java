package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class BreakpointSummaryTest {
	@Test
	public void BreakpointSummary_explicit_constructor() {
		BreakpointSummary loc = new BreakpointSummary(1, BreakendDirection.Forward, 2, 3, 4, BreakendDirection.Backward, 5, 6, 7);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
		assertEquals(7, loc.qual, 0);
	}
	@Test
	public void BreakpointSummary_BreakendSummary_constructor() {
		BreakpointSummary loc = new BreakpointSummary(
				new BreakendSummary(1, BreakendDirection.Forward, 2, 3, 8),
				new BreakendSummary(4, BreakendDirection.Backward, 5, 6, 9),
				7);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
		assertEquals(7, loc.qual, 0);
	}
}
