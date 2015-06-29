package au.edu.wehi.idsv.util;

import static org.junit.Assert.*;

import org.junit.Test;


public class IntervalUtilTest {
	private void testOverlapsWidthClosed(int expected, int s1, int e1, int s2, int e2) {
		assertEquals(expected, IntervalUtil.overlapsWidthClosed(s1, e1, s2, e2));
		assertEquals(expected, IntervalUtil.overlapsWidthClosed(s2, e2, s1, e1));
	}
	@Test
	public void overlapsWidthClosed() {
		testOverlapsWidthClosed(0, 1, 1, 2, 2);
		testOverlapsWidthClosed(0, -10, 1, 2, 10);
		testOverlapsWidthClosed(1, 1, 1, 1, 1);
		testOverlapsWidthClosed(1, 1, 1, 1, 2);
		testOverlapsWidthClosed(1, 1, 1, 0, 2);
		testOverlapsWidthClosed(2, 1, 2, 0, 2);
		testOverlapsWidthClosed(2, 1, 3, 0, 2);
		testOverlapsWidthClosed(3, 1, 3, 1, 3);
		testOverlapsWidthClosed(3, 1, 3, 0, 4);
		testOverlapsWidthClosed(3, 1, 3, 0, 4);
	}
}
