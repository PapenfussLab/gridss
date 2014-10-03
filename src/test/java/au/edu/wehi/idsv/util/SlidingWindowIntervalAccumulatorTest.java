package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class SlidingWindowIntervalAccumulatorTest {
	public SlidingWindowIntervalAccumulator accumulator;
	@Before
	public void setUp() throws Exception {
		accumulator = new SlidingWindowIntervalAccumulator();
	}
	@Test(expected=IllegalArgumentException.class)
	public void addInterval_should_require_nonnegative_interval_length() {
		accumulator.addInterval(2, 1);
	}
	@Test
	public void addInterval_should_be_start_and_end_inclusive() {
		accumulator.addInterval(1, 3);
		assertEquals(0, accumulator.getIntervalCount(0));
		assertEquals(1, accumulator.getIntervalCount(1));
		assertEquals(1, accumulator.getIntervalCount(2));
		assertEquals(1, accumulator.getIntervalCount(3));
		assertEquals(0, accumulator.getIntervalCount(4));
	}
	@Test
	public void addInterval_should_accumulate_multiple_intervals() {
		accumulator.addInterval(1, 3);
		accumulator.addInterval(2, 2);
		assertEquals(0, accumulator.getIntervalCount(0));
		assertEquals(1, accumulator.getIntervalCount(1));
		assertEquals(2, accumulator.getIntervalCount(2));
		assertEquals(1, accumulator.getIntervalCount(3));
		assertEquals(0, accumulator.getIntervalCount(4));
	}
	@Test
	public void addInterval_should_accumulate_overlapping_intervals() {
		accumulator.addInterval(1, 3);
		accumulator.addInterval(2, 4);
		assertEquals(0, accumulator.getIntervalCount(0));
		assertEquals(1, accumulator.getIntervalCount(1));
		assertEquals(2, accumulator.getIntervalCount(2));
		assertEquals(2, accumulator.getIntervalCount(3));
		assertEquals(1, accumulator.getIntervalCount(4));
		assertEquals(0, accumulator.getIntervalCount(5));
	}
	@Test(expected=IllegalArgumentException.class)
	public void addInterval_should_throw_exception_when_flushed_result_queried() {
		accumulator.addInterval(1, 3);
		accumulator.advanceTo(2);
		accumulator.getIntervalCount(1);
	}
	@Test(expected=IllegalArgumentException.class)
	public void addInterval_should_not_allow_interval_to_be_added_once_queried() {
		accumulator.getIntervalCount(1);
		accumulator.addInterval(1, 1);
	}
	@Test(expected=IllegalArgumentException.class)
	public void addInterval_should_not_allow_interval_to_be_added_once_subsequent_position_queried() {
		accumulator.getIntervalCount(5);
		accumulator.addInterval(1, 3);
	}
	@Test
	public void addInterval_should_not_interval_to_be_added_before_max_query_position() {
		accumulator.getIntervalCount(5);
		accumulator.getIntervalCount(1);
		accumulator.addInterval(6, 7);
	}
}
