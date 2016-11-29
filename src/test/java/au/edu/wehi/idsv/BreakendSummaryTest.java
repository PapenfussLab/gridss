package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class BreakendSummaryTest extends TestHelper {
	@Test
	public void BreakendSummary_Constructor() {
		BreakendSummary loc = new BreakendSummary(1, BreakendDirection.Forward, 3, 2, 4);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(4, loc.end);
		assertEquals(3, loc.nominal);
	}
	@Test
	public void containedBy_should_require_full_interval_containment() {
		assertTrue(new BreakendSummary(0, FWD, 2, 2, 3).containedBy(new BreakendSummary(0, FWD, 1, 1, 4)));
		assertTrue(new BreakendSummary(0, FWD, 2, 2, 3).containedBy(new BreakendSummary(0, FWD, 3, 2, 3)));
		assertFalse(new BreakendSummary(0, FWD, 2, 2, 3).containedBy(new BreakendSummary(0, FWD, 3, 3, 4)));
		assertFalse(new BreakendSummary(0, FWD, 2, 2, 3).containedBy(new BreakendSummary(0, FWD, 2, 1, 2)));
	}
	@Test
	public void containedBy_should_require_matching_contigs() {
		assertFalse(new BreakendSummary(0, FWD, 1, 1, 2).containedBy(new BreakendSummary(1, FWD, 1, 1, 2)));
	}
	@Test
	public void containedBy_should_require_matching_direction() {
		assertFalse(new BreakendSummary(0, FWD, 1, 1, 2).containedBy(new BreakendSummary(0, BWD, 1, 1, 2)));
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
		testOverlap(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, BWD, 1, 1, 2), null);
	}
	@Test
	public void overlaps_should_require_matching_contigs() {
		testOverlap(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(1, FWD, 1, 1, 2), null);
	}
	@Test
	public void overlaps_should_require_at_least_one_position_in_common() {
		testOverlap(new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 1, 1, 4), new BreakendSummary(0, FWD, 2, 2, 3));
		testOverlap(new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 2, 2, 3));
		testOverlap(new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 3, 3, 4), new BreakendSummary(0, FWD, 3, 3, 3));
		testOverlap(new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, FWD, 2, 2, 2));
		testOverlap(new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 1, 1, 1), null);
		testOverlap(new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 4, 4, 4), null);
	}
	@Test
	public void equals_should_require_full_match() {
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, FWD, 1, 1, 2));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(1, FWD, 1, 1, 2));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, BWD, 2, 2, 2));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, FWD, 1, 1, 3));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, FWD, 1, 1, 3));
		assertNotEquals(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, FWD, 2, 1, 2));
	}
	@Test
	public void expandBounds_should_extend_both_sides() {
		assertEquals(new BreakendSummary(0, FWD, 4, 3, 8), new BreakendSummary(0, FWD, 4, 4, 7).expandBounds(1));
		assertEquals(new BreakendSummary(0, FWD, 4, 2, 9), new BreakendSummary(0, FWD, 4, 4, 7).expandBounds(2));
	}
	@Test
	public void expandBounds_should_extend_past_chr_bounds() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("test", 10));
		assertEquals(new BreakendSummary(0, FWD, 9, 7, 11), new BreakendSummary(0, FWD, 9, 9, 9).expandBounds(2));
		assertEquals(new BreakendSummary(0, FWD, 1, -1, 3), new BreakendSummary(0, FWD, 1, 1, 1).expandBounds(2));
	}
	@Test
	public void compressBounds_should_reduce_both_sides() {
		assertEquals(new BreakendSummary(0, FWD, 4, 4, 6), new BreakendSummary(0, FWD, 4, 3, 7).compressBounds(1));
		assertEquals(new BreakendSummary(0, FWD, 5, 5, 5), new BreakendSummary(0, FWD, 5, 3, 7).compressBounds(2));
	}
	@Test
	public void compressBounds_should_not_reduce_to_less_than_single_base() {
		assertEquals(new BreakendSummary(0, FWD, 5, 5, 5), new BreakendSummary(0, FWD, 5, 3, 7).compressBounds(10));
	}
	@Test
	public void compressBounds_should_move_nominal_to_within_valid_range() {
		assertEquals(new BreakendSummary(0, FWD, 5, 5, 5), new BreakendSummary(0, FWD, 3, 3, 7).compressBounds(10));
	}
	@Test
	public void compressBounds_should_reduce_to_lower_middle_base() {
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 1), new BreakendSummary(0, FWD, 1, 1, 1).compressBounds(0));
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 2), new BreakendSummary(0, FWD, 1, 1, 2).compressBounds(0));
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 3), new BreakendSummary(0, FWD, 1, 1, 3).compressBounds(0));
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 4), new BreakendSummary(0, FWD, 1, 1, 4).compressBounds(0));
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 5), new BreakendSummary(0, FWD, 1, 1, 5).compressBounds(0));
		
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 1), new BreakendSummary(0, FWD, 1, 1, 1).compressBounds(1));
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 1), new BreakendSummary(0, FWD, 1, 1, 2).compressBounds(1));
		assertEquals(new BreakendSummary(0, FWD, 2, 2, 2), new BreakendSummary(0, FWD, 2, 1, 3).compressBounds(1));
		assertEquals(new BreakendSummary(0, FWD, 2, 2, 3), new BreakendSummary(0, FWD, 2, 1, 4).compressBounds(1));
		assertEquals(new BreakendSummary(0, FWD, 2, 2, 4), new BreakendSummary(0, FWD, 2, 1, 5).compressBounds(1));
		
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 1), new BreakendSummary(0, FWD, 1, 1, 1).compressBounds(2));
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 1), new BreakendSummary(0, FWD, 1, 1, 2).compressBounds(2));
		assertEquals(new BreakendSummary(0, FWD, 2, 2, 2), new BreakendSummary(0, FWD, 2, 1, 3).compressBounds(2));
		assertEquals(new BreakendSummary(0, FWD, 2, 2, 2), new BreakendSummary(0, FWD, 2, 1, 4).compressBounds(2));
		assertEquals(new BreakendSummary(0, FWD, 3, 3, 3), new BreakendSummary(0, FWD, 3, 1, 5).compressBounds(2));
		assertEquals(new BreakendSummary(0, FWD, 3, 3, 4), new BreakendSummary(0, FWD, 3, 1, 6).compressBounds(2));
		assertEquals(new BreakendSummary(0, FWD, 3, 3, 5), new BreakendSummary(0, FWD, 3, 1, 7).compressBounds(2));
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_ensure_referenceindex_bounds() {
		new BreakendSummary(-1, FWD, 1, 1, 1);
	}
	@Test
	public void constructor_should_not_ensure_position_bounds() {
		new BreakendSummary(0, FWD, 0, 0, 1);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_ensure_positive_interval() {
		new BreakendSummary(0, FWD, 1, 2, 1);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_ensure_nominal_position_is_valid() {
		new BreakendSummary(0, FWD, 3, 1, 2);
	}
	@Test
	public void isValid_should_check_bounds() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("test", 10));
		dict.addSequence(new SAMSequenceRecord("test2", 20));
		assertTrue(new BreakendSummary(0, FWD, 1, 1, 1).isValid(dict));
		assertTrue(new BreakendSummary(0, FWD, 10, 10, 10).isValid(dict));
		assertTrue(new BreakendSummary(0, FWD, 1, 1, 10).isValid(dict));
		assertTrue(new BreakendSummary(1, FWD, 1, 1, 20).isValid(dict));
		
		assertFalse(new BreakendSummary(2, FWD, 1, 1, 1).isValid(dict));
		assertFalse(new BreakendSummary(0, FWD, 1, 1, 11).isValid(dict));
	}
	@Test
	public void isValid_should_check_start_end() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("test", 10));
		dict.addSequence(new SAMSequenceRecord("test2", 10));
		assertTrue(new BreakendSummary(0, FWD, 2, 2, 2).isValid(dict));
	}
	@Test
	public void getNominalPosition_should_report_nominal_position_only() {
		assertEquals(new BreakendSummary(0, FWD, 15), new BreakendSummary(0, FWD, 15, 1, 100).getNominalPosition());
	}
	@Test
	public void asValidFor_should_adjust_bounds() {
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		dict.addSequence(new SAMSequenceRecord("test", 10));
		assertEquals(new BreakendSummary(0, FWD, 4, 1, 10), new BreakendSummary(0, FWD, 4, -10, 20).asValidFor(dict));
		assertEquals(new BreakendSummary(0, FWD, 10, 10, 10), new BreakendSummary(0, FWD, 16, 15, 20).asValidFor(dict));
		assertEquals(new BreakendSummary(0, FWD, 1, 1, 1), new BreakendSummary(0, FWD, -10, -10, -10).asValidFor(dict));
	}
	@Test
	public void adjustPosition_should_move_according_to_breakend_direction() {
		assertEquals(new BreakendSummary(0, FWD, 2, 0, 5), new BreakendSummary(0, FWD, 2, 1, 3).adjustPosition(1, 2, false));
		assertEquals(new BreakendSummary(0, BWD, 2, -1, 4), new BreakendSummary(0, BWD, 2, 1, 3).adjustPosition(1, 2, false));
	}
	@Test
	public void adjustPosition_should_move_nominal_by_adjustment() {
		assertEquals(new BreakendSummary(0, FWD, 2, -1, 5), new BreakendSummary(0, FWD, 2, 1, 3).adjustPosition(2, 2, true));
		assertEquals(new BreakendSummary(0, BWD, 3, 1, 5), new BreakendSummary(0, BWD, 2, 1, 3).adjustPosition(2, 0, true));
	}
	@Test
	public void adjustPosition_should_retain_nominal_position() {
		assertEquals(new BreakendSummary(0, FWD, 10, 10, 20), new BreakendSummary(0, FWD, 10).adjustPosition(0, 10, false));
	}
}
