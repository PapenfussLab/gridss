package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class MathUtilTest {
	@Test
	public void phredOr_should_be_exact_for_single_score() {
		assertEquals(5, MathUtil.phredOr(5), 0);
	}
	@Test
	public void phredOr_should_minimise_error() {
		// http://www.wolframalpha.com/input/?i=-10+*+%28log%2810%2C+1+-+%281+-+10%5E%28-5%2F10%29%29+*+%281+-+10%5E%28-5%2F10%29%29%29+%29
		assertEquals(2.737166564315747414571566754814524996985593328024647081266066, MathUtil.phredOr(5, 5), 0);
		assertEquals(4.999906094489629897297010515927232501782840270893945897373360, MathUtil.phredOr(5, 50), 0);
		assertEquals(0.999998875501369913289357102130159464124571763746365284402279, MathUtil.phredOr(1, 60), 0);
		assertEquals(0.99999999988755012243472182763604174153706462136143019889429, MathUtil.phredOr(1, 100), 0.000000000000001);
	}
}
