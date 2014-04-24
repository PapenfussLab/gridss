package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class BreakpointIntervalTest {
	@Test
	public void BreakpointInterval_explicit_constructor() {
		BreakpointInterval loc = new BreakpointInterval(1, BreakpointDirection.Forward, 2, 3, 4, BreakpointDirection.Backward, 5, 6, 7);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakpointDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakpointDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
		assertEquals(7, loc.qual, 0);
	}
	@Test
	public void BreakpointInterval_BreakpointLocation_constructor() {
		BreakpointInterval loc = new BreakpointInterval(
				new BreakpointLocation(1, BreakpointDirection.Forward, 2, 3, 8),
				new BreakpointLocation(4, BreakpointDirection.Backward, 5, 6, 9),
				7);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakpointDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakpointDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
		assertEquals(7, loc.qual, 0);
	}
}
