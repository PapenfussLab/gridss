package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public class ArrayHelperTest {
	@Test
	public void should_add_arrays() {
		assertArrayEquals(new int[] { 4, 6 }, ArrayHelper.add(new int[] { 1, 2 }, new int[] { 3, 4 }));
		assertArrayEquals(new int[] { 3, 4 }, ArrayHelper.add(null, new int[] { 3, 4 }));
		assertArrayEquals(new int[] { 3, 4 }, ArrayHelper.add(new int[] { 3, 4 }, null));
		assertArrayEquals(new int[] { }, ArrayHelper.add(null, null));
		assertArrayEquals(new int[] { 4, 6, 5 }, ArrayHelper.add(new int[] { 1, 2 }, new int[] { 3, 4, 5 }));
		assertArrayEquals(new int[] { 5, 7, 3 }, ArrayHelper.add(new int[] { 1, 2, 3 }, new int[] { 4, 5 }));
	}
	@Test
	public void should_subtract_arrays() {
		assertArrayEquals(new int[] { -2, -2 }, ArrayHelper.subtract(new int[] { 1, 2 }, new int[] { 3, 4 }));
		assertArrayEquals(new int[] { -3, -4 }, ArrayHelper.subtract(null, new int[] { 3, 4 }));
		assertArrayEquals(new int[] { 3, 4 }, ArrayHelper.subtract(new int[] { 3, 4 }, null));
		assertArrayEquals(new int[] { }, ArrayHelper.subtract(null, null));
		assertArrayEquals(new int[] { -2, -2, -5 }, ArrayHelper.subtract(new int[] { 1, 2 }, new int[] { 3, 4, 5 }));
		assertArrayEquals(new int[] { -3, -3, 3 }, ArrayHelper.subtract(new int[] { 1, 2, 3 }, new int[] { 4, 5 }));
	}
	@Test
	public void should_add_position_value() {
		assertArrayEquals(new int[] { 1 }, ArrayHelper.add(null, 0, 1));
		assertArrayEquals(new int[] { 0, 2 }, ArrayHelper.add(null, 1, 2));
		assertArrayEquals(new int[] { 0, 2 }, ArrayHelper.add(new int[] { }, 1, 2));
		assertArrayEquals(new int[] { 0, 2 }, ArrayHelper.add(new int[] {0 }, 1, 2));
		assertArrayEquals(new int[] { 0, 2, 0 }, ArrayHelper.add(new int[] {0, 0, 0 }, 1, 2));
		assertArrayEquals(new int[] { 4, 7, 6 }, ArrayHelper.add(new int[] {4, 5, 6 }, 1, 2));
	}
}
