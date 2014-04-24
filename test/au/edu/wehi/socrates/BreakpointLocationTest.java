package au.edu.wehi.socrates;

import static org.junit.Assert.*;

import org.junit.Test;

public class BreakpointLocationTest {
	@Test
	public void BreakpointLocation_Constructor() {
		BreakpointLocation loc = new BreakpointLocation(1, BreakpointDirection.Forward, 2, 3, 4);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakpointDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.qual, 0);
	}
}
