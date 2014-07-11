package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class BreakendSummaryTest extends TestHelper {
	@Test
	public void BreakendSummary_Constructor() {
		EvidenceMetrics em = new EvidenceMetrics();
		BreakendSummary loc = new BreakendSummary(1, BreakendDirection.Forward, 2, 3, em);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(em, loc.evidence);
	}
	@Test
	public void containedBy_should_require_full_interval_containment() {
		assertTrue(new BreakendSummary(0, FWD, 2, 3, null).containedBy(new BreakendSummary(0, FWD, 1, 4, null)));
		assertTrue(new BreakendSummary(0, FWD, 2, 3, null).containedBy(new BreakendSummary(0, FWD, 2, 3, null)));
		assertFalse(new BreakendSummary(0, FWD, 2, 3, null).containedBy(new BreakendSummary(0, FWD, 3, 4, null)));
		assertFalse(new BreakendSummary(0, FWD, 2, 3, null).containedBy(new BreakendSummary(0, FWD, 1, 2, null)));
	}
	@Test
	public void containedBy_should_require_matching_contigs() {
		assertFalse(new BreakendSummary(0, FWD, 1, 2, null).containedBy(new BreakendSummary(1, FWD, 1, 2, null)));
	}
	@Test
	public void containedBy_should_require_matching_direction() {
		assertFalse(new BreakendSummary(0, FWD, 1, 2, null).containedBy(new BreakendSummary(0, BWD, 1, 2, null)));
	}
	/**
	 * Overlaps should be Commutative
	 * @param a
	 * @param b
	 */
	public void testOverlap(BreakendSummary a, BreakendSummary b, boolean expected) {
		assertEquals(a.overlaps(b), expected);
		assertEquals(b.overlaps(a), expected);
	}
	@Test
	public void overlaps_should_require_matching_direction() {
		testOverlap(new BreakendSummary(0, FWD, 1, 2, null), new BreakendSummary(0, BWD, 1, 2, null), false);
	}
	@Test
	public void overlaps_should_require_matching_contigs() {
		testOverlap(new BreakendSummary(0, FWD, 1, 2, null), new BreakendSummary(1, FWD, 1, 2, null), false);
	}
	@Test
	public void overlaps_should_require_at_least_one_position_in_common() {
		testOverlap(new BreakendSummary(0, FWD, 2, 3, null), new BreakendSummary(0, FWD, 1, 4, null), true);
		testOverlap(new BreakendSummary(0, FWD, 2, 3, null), new BreakendSummary(0, FWD, 2, 3, null), true);
		testOverlap(new BreakendSummary(0, FWD, 2, 3, null), new BreakendSummary(0, FWD, 3, 4, null), true);
		testOverlap(new BreakendSummary(0, FWD, 2, 3, null), new BreakendSummary(0, FWD, 1, 2, null), true);
		testOverlap(new BreakendSummary(0, FWD, 2, 3, null), new BreakendSummary(0, FWD, 1, 1, null), false);
		testOverlap(new BreakendSummary(0, FWD, 2, 3, null), new BreakendSummary(0, FWD, 4, 4, null), false);
	}
}
