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
	public void overlaps_should_consider_breakpoint_remote_location_when_downcast() {
		BreakendSummary be = new BreakpointSummary(0, FWD, 2, 3, 1, FWD, 2, 3);
		assertTrue(be.overlaps(new BreakendSummary(0, FWD, 3, 3)));
		assertFalse(be.overlaps(new BreakpointSummary(0, FWD, 3, 3, 0, FWD, 1, 1)));
		assertTrue(be.overlaps(new BreakpointSummary(0, FWD, 3, 3, 1, FWD, 2, 2)));
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
		assertEquals(new BreakpointSummary(0, FWD, 3, 8, 1, BWD, 9, 12), new BreakpointSummary(0, FWD, 4, 7, 1, BWD, 10, 11).expandBounds(1));
	}
	@Test
	public void compressBounds_should_compress_both_sides() {
		assertEquals(new BreakpointSummary(0, FWD, 5, 6, 1, BWD, 11, 11), new BreakpointSummary(0, FWD, 4, 7, 1, BWD, 10, 12).compressBounds(1));
		assertEquals(new BreakpointSummary(0, FWD, 5, 6, 1, BWD, 11, 11), ((BreakendSummary)new BreakpointSummary(0, FWD, 4, 7, 1, BWD, 10, 12)).compressBounds(1));
	}
	@Test
	public void compressBounds_should_reduce_lower_coordinate_down_upper_coordinate_up() {
		assertEquals(new BreakpointSummary(0, FWD, 4, 4, 1, FWD, 5, 5), new BreakpointSummary(0, FWD, 4, 5, 1, FWD, 4, 5).compressBounds(1));
		assertEquals(new BreakpointSummary(1, FWD, 5, 5, 0, FWD, 4, 4), new BreakpointSummary(1, FWD, 4, 5, 0, FWD, 4, 5).compressBounds(1));
		assertEquals(new BreakpointSummary(0, FWD, 4, 4, 0, FWD, 6, 6), new BreakpointSummary(0, FWD, 4, 5, 0, FWD, 5, 6).compressBounds(5));
		assertEquals(new BreakpointSummary(0, FWD, 6, 6, 0, FWD, 4, 4), new BreakpointSummary(0, FWD, 5, 6, 0, FWD, 4, 5).compressBounds(5));
		assertEquals(new BreakpointSummary(0, FWD, 3, 3, 0, FWD, 1, 1), new BreakpointSummary(0, FWD, 1, 4, 0, FWD, 1, 2).compressBounds(5));
	}
	@Test
	public void compressBounds_should_reduce_both_down_if_symmetrical() {
		assertEquals(new BreakpointSummary(0, FWD, 1, 1, 0, FWD, 1, 1), new BreakpointSummary(0, FWD, 1, 2, 0, FWD, 1, 2).compressBounds(1));
	}
	@Test
	public void isHighBreakend() {
		assertFalse(new BreakpointSummary(1, FWD, 10, 20, 1, FWD, 10, 20).isHighBreakend());
		assertFalse(new BreakpointSummary(1, FWD, 10, 20, 1, BWD, 10, 20).isHighBreakend());
		assertFalse(new BreakpointSummary(0, FWD, 10, 20, 1, FWD, 10, 20).isHighBreakend());
		assertFalse(new BreakpointSummary(1, FWD, 9, 20, 1, FWD, 10, 20).isHighBreakend());
		assertFalse(new BreakpointSummary(1, FWD, 10, 19, 1, FWD, 10, 20).isHighBreakend());
		assertTrue(new BreakpointSummary(2, FWD, 10, 20, 1, FWD, 10, 20).isHighBreakend());
		assertTrue(new BreakpointSummary(1, FWD, 11, 20, 1, FWD, 10, 20).isHighBreakend());
		assertTrue(new BreakpointSummary(1, FWD, 10, 21, 1, FWD, 10, 20).isHighBreakend());
	}
	@Test
	public void isLowBreakend() {
		assertFalse(new BreakpointSummary(1, FWD, 10, 20, 1, FWD, 10, 20).isLowBreakend());
		assertFalse(new BreakpointSummary(1, FWD, 10, 20, 1, BWD, 10, 20).isLowBreakend());
		assertTrue(new BreakpointSummary(0, FWD, 10, 20, 1, FWD, 10, 20).isLowBreakend());
		assertTrue(new BreakpointSummary(1, FWD, 9, 20, 1, FWD, 10, 20).isLowBreakend());
		assertTrue(new BreakpointSummary(1, FWD, 10, 19, 1, FWD, 10, 20).isLowBreakend());
		assertFalse(new BreakpointSummary(2, FWD, 10, 20, 1, FWD, 10, 20).isLowBreakend());
		assertFalse(new BreakpointSummary(1, FWD, 11, 20, 1, FWD, 10, 20).isLowBreakend());
		assertFalse(new BreakpointSummary(1, FWD, 10, 21, 1, FWD, 10, 20).isLowBreakend());
	}
	@Test
	public void couldBeDeletionOfSize() {
		assertTrue(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 3, 3).couldBeDeletionOfSize(1, 1));
		assertFalse(new BreakpointSummary(0, FWD, 1, 1, 0, FWD, 3, 3).couldBeDeletionOfSize(1, 1));
		assertFalse(new BreakpointSummary(0, BWD, 1, 1, 0, BWD, 3, 3).couldBeDeletionOfSize(1, 1));
		assertFalse(new BreakpointSummary(0, BWD, 1, 1, 0, FWD, 3, 3).couldBeDeletionOfSize(1, 1));
		
		assertTrue(new BreakpointSummary(0, BWD, 3, 3, 0, FWD, 1, 1).couldBeDeletionOfSize(1, 1));
		assertTrue(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 3).couldBeDeletionOfSize(1, 1));
		assertTrue(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 4, 4).couldBeDeletionOfSize(2, 2));
		assertTrue(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 3, 3).couldBeDeletionOfSize(1, 2));
		assertTrue(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 4, 4).couldBeDeletionOfSize(1, 2));
	}
	@Test
	public void remoteBreakpoint() {
		assertEquals(new BreakpointSummary(4, BWD, 5, 6, 1, FWD, 2, 3), new BreakpointSummary(1, FWD, 2, 3, 4, BWD, 5, 6).remoteBreakpoint());
	}
	@Test
	public void getEventSize_deletion() {
		assertEquals(0, (int)new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 2, 2).getEventSize());
		assertEquals(0, (int)new BreakpointSummary(0, BWD, 2, 2, 0, FWD, 1, 1).getEventSize());
		assertEquals(1, (int)new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 3, 3).getEventSize());
		assertEquals(1, (int)new BreakpointSummary(0, BWD, 3, 3, 0, FWD, 1, 1).getEventSize());
	}
	@Test
	public void getEventSize_duplication() {
		assertEquals(1, (int)new BreakpointSummary(0, FWD, 2, 2, 0, BWD, 2, 2).getEventSize());
		assertEquals(2, (int)new BreakpointSummary(0, FWD, 3, 3, 0, BWD, 2, 2).getEventSize());
		assertEquals(1, (int)new BreakpointSummary(0, BWD, 2, 2, 0, FWD, 2, 2).getEventSize());
		assertEquals(2, (int)new BreakpointSummary(0, BWD, 2, 2, 0, FWD, 3, 3).getEventSize());
	}
	@Test
	public void getEventSize_inversion() {
		assertEquals(2, (int)new BreakpointSummary(0, FWD, 2, 2, 0, FWD, 4, 4).getEventSize());
		assertEquals(2, (int)new BreakpointSummary(0, BWD, 3, 3, 0, BWD, 5, 5).getEventSize());
		assertEquals(2, (int)new BreakpointSummary(0, FWD, 4, 4, 0, FWD, 2, 2).getEventSize());
		assertEquals(2, (int)new BreakpointSummary(0, BWD, 3, 3, 0, BWD, 5, 5).getEventSize());
	}
}
