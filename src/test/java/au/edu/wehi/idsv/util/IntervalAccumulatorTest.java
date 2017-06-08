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
	public void bins_should_span_interval() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 3, 1);
		ia.finaliseBins();
		int[] binStarts = ia.getBinStarts().toIntArray();
		Assert.assertArrayEquals(new int[] { 1, 2, 3}, binStarts);
		Assert.assertArrayEquals(new int[] { 1, 1, 1}, IntStream.of(ia.getBinStarts().toIntArray()).map(i -> ia.getBinSize(i)).toArray());
	}
	@Test
	public void should_average_across_bin() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 5, 10);
		ia.finaliseBins();
		ia.add(1, 2, 1);
		Assert.assertEquals(0.4, ia.getMeanValue(1), 0);
	}
	@Test
	public void should_average_across_adjacent_bins() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 5, 2);
		ia.finaliseBins();
		ia.add(2, 3, 1);
		Assert.assertEquals(0.5, ia.getMeanValue(1), 0);
		Assert.assertEquals(0.5, ia.getMeanValue(3), 0);
	}
	@Test
	public void should_average_across_multiple_bins() {
		IntervalAccumulator ia = new IntervalAccumulator(1, 5, 2);
		ia.finaliseBins();
		// 12345
		// 11223 bin
		//  **** coverage
		//  **
		ia.add(2, 5, 1);
		ia.add(2, 3, 1);
		Assert.assertEquals(1, ia.getMeanValue(1), 0);
		Assert.assertEquals(1.5, ia.getMeanValue(3), 0);
		Assert.assertEquals(1, ia.getMeanValue(5), 0);
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
