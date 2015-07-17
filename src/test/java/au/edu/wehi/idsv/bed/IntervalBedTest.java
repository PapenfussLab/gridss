package au.edu.wehi.idsv.bed;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class IntervalBedTest extends TestHelper {
	@Test
	public void should_round_trip() {
		IntervalBed bed = new IntervalBed(getContext().getDictionary(), getContext().getLinear());
		bed.addInterval(1, 3, 5);
	}
	@Test
	public void overlap_should_return_true_if_at_least_one_base_overlaps() {
		IntervalBed bed = new IntervalBed(getContext().getDictionary(), getContext().getLinear());
		bed.addInterval(1, 3, 5);
		bed.addInterval(1, 7, 9);
		assertFalse(bed.overlaps(0, 3, 5));
		assertFalse(bed.overlaps(1, 2, 2));
		assertFalse(bed.overlaps(1, 6, 6));
		assertTrue(bed.overlaps(1, 3, 5));
		assertTrue(bed.overlaps(1, 3, 3));
		assertTrue(bed.overlaps(1, 3, 4));
		assertTrue(bed.overlaps(1, 4, 4));
		assertTrue(bed.overlaps(1, 5, 5));
		assertTrue(bed.overlaps(1, 5, 6));
		assertTrue(bed.overlaps(1, 2, 3));
		assertTrue(bed.overlaps(1, 2, 4));
		assertTrue(bed.overlaps(1, 2, 5));
		assertTrue(bed.overlaps(1, 2, 10));
		assertTrue(bed.overlaps(1, 1, 11));
	}
}
