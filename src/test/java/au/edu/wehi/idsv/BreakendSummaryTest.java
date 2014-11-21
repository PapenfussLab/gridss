package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class BreakendSummaryTest extends TestHelper {
	@Test
	public void BreakendSummary_Constructor() {
		BreakendSummary loc = new BreakendSummary(1, BreakendDirection.Forward, 2, 3);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
	}
	@Test
	public void containedBy_should_require_full_interval_containment() {
		assertTrue(new BreakendSummary(0, FWD, 2, 3).containedBy(new BreakendSummary(0, FWD, 1, 4)));
		assertTrue(new BreakendSummary(0, FWD, 2, 3).containedBy(new BreakendSummary(0, FWD, 2, 3)));
		assertFalse(new BreakendSummary(0, FWD, 2, 3).containedBy(new BreakendSummary(0, FWD, 3, 4)));
		assertFalse(new BreakendSummary(0, FWD, 2, 3).containedBy(new BreakendSummary(0, FWD, 1, 2)));
	}
	@Test
	public void containedBy_should_require_matching_contigs() {
		assertFalse(new BreakendSummary(0, FWD, 1, 2).containedBy(new BreakendSummary(1, FWD, 1, 2)));
	}
	@Test
	public void containedBy_should_require_matching_direction() {
		assertFalse(new BreakendSummary(0, FWD, 1, 2).containedBy(new BreakendSummary(0, BWD, 1, 2)));
	}
	/**
	 * Overlaps should be Commutative
	 * @param a
	 * @param b
	 */
	public void testOverlap(BreakendSummary a, BreakendSummary b, BreakendSummary expected) {
		assertEquals(a.overlaps(b), expected != null);
		assertEquals(b.overlaps(a), expected != null);
		assertEquals(BreakendSummary.overlapOf(a, b), expected);
		assertEquals(BreakendSummary.overlapOf(b, a), expected);
	}
	@Test
	public void overlaps_should_require_matching_direction() {
		testOverlap(new BreakendSummary(0, FWD, 1, 2), new BreakendSummary(0, BWD, 1, 2), null);
	}
	@Test
	public void overlaps_should_require_matching_contigs() {
		testOverlap(new BreakendSummary(0, FWD, 1, 2), new BreakendSummary(1, FWD, 1, 2), null);
	}
	@Test
	public void overlaps_should_require_at_least_one_position_in_common() {
		testOverlap(new BreakendSummary(0, FWD, 2, 3), new BreakendSummary(0, FWD, 1, 4), new BreakendSummary(0, FWD, 2, 3));
		testOverlap(new BreakendSummary(0, FWD, 2, 3), new BreakendSummary(0, FWD, 2, 3), new BreakendSummary(0, FWD, 2, 3));
		testOverlap(new BreakendSummary(0, FWD, 2, 3), new BreakendSummary(0, FWD, 3, 4), new BreakendSummary(0, FWD, 3, 3));
		testOverlap(new BreakendSummary(0, FWD, 2, 3), new BreakendSummary(0, FWD, 1, 2), new BreakendSummary(0, FWD, 2, 2));
		testOverlap(new BreakendSummary(0, FWD, 2, 3), new BreakendSummary(0, FWD, 1, 1), null);
		testOverlap(new BreakendSummary(0, FWD, 2, 3), new BreakendSummary(0, FWD, 4, 4), null);
	}
	@Test
	public void equals_should_require_full_match() {
		assertEquals(new BreakendSummary(0, FWD, 1, 2), new BreakendSummary(0, FWD, 1, 2));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 2), new BreakendSummary(1, FWD, 1, 2));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 2), new BreakendSummary(0, BWD, 2, 2));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 2), new BreakendSummary(0, FWD, 1, 3));
	}
}
