package au.edu.wehi.idsv.util;

import java.util.stream.IntStream;

import org.junit.Assert;
import org.junit.Test;

public class IntervalAccumulatorTest {
	@Test
	public void should_split_into_bins() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 5, 2);
		ia.finaliseBins();
		int[] binStarts = ia.getBinStarts().toIntArray();
		Assert.assertArrayEquals(new int[] { 1, 3, 5}, binStarts);
		Assert.assertArrayEquals(new int[] { 2, 2, 1}, IntStream.of(ia.getBinStarts().toIntArray()).map(i -> ia.getBinSize(i)).toArray());
	}
	@Test
	public void should_pro_rata_across_adjacent_bins() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 5, 2);
		ia.finaliseBins();
		ia.add(2, 3, 1);
		Assert.assertEquals(0.5, ia.getValue(1), 0);
		Assert.assertEquals(0.5, ia.getValue(3), 0);
	}
	@Test
	public void should_pro_rata_across_multiple_bins() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 5, 2);
		ia.finaliseBins();
		ia.add(2, 5, 1);
		Assert.assertEquals(0.25, ia.getValue(1), 0);
		Assert.assertEquals(0.5, ia.getValue(3), 0);
		Assert.assertEquals(0.25, ia.getValue(5), 0);
	}
	@Test
	public void should_split_bins_at_position() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 5, 3);
		ia.splitBin(2);
		ia.finaliseBins();
		Assert.assertArrayEquals(new int[] { 1, 2, 4}, ia.getBinStarts().toIntArray());
		Assert.assertArrayEquals(new int[] { 1, 2, 2}, IntStream.of(ia.getBinStarts().toIntArray()).map(i -> ia.getBinSize(i)).toArray());
	}
}
