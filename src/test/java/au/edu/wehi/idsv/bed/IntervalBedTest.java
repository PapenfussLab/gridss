package au.edu.wehi.idsv.bed;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.QueryInterval;


public class IntervalBedTest extends TestHelper {
	@Test
	public void should_round_trip() throws IOException {
		IntervalBed bed = new IntervalBed(getContext().getLinear());
		bed.addInterval(1, 3, 5);
		TemporaryFolder folder = new TemporaryFolder();
		folder.create();
		File f = folder.newFile("IntervalBedTest.bed");
		f.delete();
		bed.write(f, "IntervalBedTest");
		IntervalBed bed2 = new IntervalBed(getContext().getLinear(), f);
		assertFalse(bed2.overlaps(1, 2, 2));
		assertTrue(bed2.overlaps(1, 3, 3));
		assertTrue(bed2.overlaps(1, 4, 4));
		assertTrue(bed2.overlaps(1, 5, 5));
		assertFalse(bed2.overlaps(1, 6, 6));
		f.delete();
		folder.delete();
	}
	@Test
	public void overlap_should_return_true_if_at_least_one_base_overlaps() {
		IntervalBed bed = new IntervalBed(getContext().getLinear());
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
	@Test
	public void remove_should_split_interval() {
		IntervalBed bed = new IntervalBed(getContext().getLinear());
		bed.addInterval(0, 1, 10);
		IntervalBed toRemove = new IntervalBed(getContext().getLinear());
		toRemove.addInterval(0, 1, 1);
		toRemove.addInterval(1, 1, 1);
		toRemove.addInterval(0, 5, 6);
		bed.remove(toRemove);
		// 1234567890
		// X   XX
		QueryInterval[] qi = bed.asQueryInterval();
		assertEquals(2, qi.length);
		assertEquals(0, qi[0].referenceIndex);
		assertEquals(2, qi[0].start);
		assertEquals(4, qi[0].end);
		assertEquals(0, qi[1].referenceIndex);
		assertEquals(7, qi[1].start);
		assertEquals(10, qi[1].end);
	}
	@Test
	public void QueryInterval_should_round_trip() {
		QueryInterval[] qi = new QueryInterval[] {
			new QueryInterval(0, 1, 2),
			new QueryInterval(2, 4, 8),
			new QueryInterval(2, 10, 10),
		};
		IntervalBed bed = new IntervalBed(getContext().getLinear(), qi);
		QueryInterval[] result = bed.asQueryInterval();
		assertEquals(qi.length, result.length);
		for (int i = 0; i < qi.length; i++) {
			assertEquals(qi[i].referenceIndex, result[i].referenceIndex);
			assertEquals(qi[i].start, result[i].start);
			assertEquals(qi[i].end, result[i].end);
		}
	}
}
