package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;


public class BreakpointSummaryTest extends TestHelper {
	@Test
	public void BreakpointSummary_explicit_constructor() {
		BreakpointSummary loc = new BreakpointSummary(1, BreakendDirection.Forward, 2, 3, 4, BreakendDirection.Backward, 5, 6);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
	}
	@Test
	public void BreakpointSummary_BreakendSummary_constructor() {
		BreakpointSummary loc = new BreakpointSummary(
				new BreakendSummary(1, BreakendDirection.Forward, 2, 3),
				new BreakendSummary(4, BreakendDirection.Backward, 5, 6));
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
	}
	@Test
	public void overlaps_should_consider_breakpoint_remote_location() {
		BreakpointSummary bp = new BreakpointSummary(0, FWD, 2, 3, 1, FWD, 2, 3); 
		for (BreakendSummary local : new BreakendSummary[] { 
				new BreakendSummary(0, FWD, 1, 4),
				new BreakendSummary(0, FWD, 2, 3),
				new BreakendSummary(0, FWD, 3, 4),
				new BreakendSummary(0, FWD, 1, 2)}) {
			for (BreakendSummary remote : new BreakendSummary[] { 
					new BreakendSummary(1, FWD, 1, 4),
					new BreakendSummary(1, FWD, 2, 3),
					new BreakendSummary(1, FWD, 3, 4),
					new BreakendSummary(1, FWD, 1, 2)}) {
				assertTrue(new BreakpointSummary(local, remote).overlaps(bp));
				assertTrue(bp.overlaps(new BreakpointSummary(local, remote)));
				// should not flip local and remote coordinates when checking for a match
				assertFalse(new BreakpointSummary(remote, local).overlaps(bp));
				assertFalse(bp.overlaps(new BreakpointSummary(remote, local)));
			}
		}
		// local does not overlap
		for (BreakendSummary local : new BreakendSummary[] { 
				new BreakendSummary(0, FWD, 1, 1),
				new BreakendSummary(0, FWD, 4, 4)}) {
			for (BreakendSummary remote : new BreakendSummary[] { 
					new BreakendSummary(1, FWD, 1, 4),
					new BreakendSummary(1, FWD, 2, 3),
					new BreakendSummary(1, FWD, 3, 4),
					new BreakendSummary(1, FWD, 1, 2),
					new BreakendSummary(1, FWD, 1, 1),
					new BreakendSummary(1, FWD, 4, 4)}) {
				assertFalse(new BreakpointSummary(local, remote).overlaps(bp));
				assertFalse(bp.overlaps(new BreakpointSummary(local, remote)));
			}
		}
		// remote does not overlap
		for (BreakendSummary local : new BreakendSummary[] { 
				new BreakendSummary(0, FWD, 1, 4),
				new BreakendSummary(0, FWD, 2, 3),
				new BreakendSummary(0, FWD, 3, 4),
				new BreakendSummary(0, FWD, 1, 2),
				new BreakendSummary(0, FWD, 1, 1),
				new BreakendSummary(0, FWD, 4, 4)}) {
			for (BreakendSummary remote : new BreakendSummary[] { 
					new BreakendSummary(1, FWD, 1, 1),
					new BreakendSummary(1, FWD, 4, 4)}) {
				assertFalse(new BreakpointSummary(local, remote).overlaps(bp));
				assertFalse(bp.overlaps(new BreakpointSummary(local, remote)));
			}
		}
	}
	@Test
	public void expandBounds_should_expand_both_sides() {
		assertEquals(new BreakpointSummary(0, FWD, 3, 8, 1, BWD, 9, 12), new BreakpointSummary(0, FWD, 4, 7, 1, BWD, 10, 11).expandBounds(1, SMALL_FA.getSequenceDictionary()));
	}
}
