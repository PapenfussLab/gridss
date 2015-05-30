package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class TestHelperTest extends TestHelper {
	@Test
	public void test_getRef() {
		assertEquals("CATT", S(getRef(RANDOM, 4, 5, FWD, 5)));
	}
	@Test
	public void test_getBasePos() {
		assertEquals(1, getBasePos(1, FWD, 1));
		
		assertEquals(10, getBasePos(10, FWD, 1));
		assertEquals(9, getBasePos(10, FWD, 2));
		assertEquals(8, getBasePos(10, FWD, 3));
		
		assertEquals(10, getBasePos(10, BWD, 1));
		assertEquals(11, getBasePos(10, BWD, 2));
		assertEquals(12, getBasePos(10, BWD, 3));
	}
	@Test
	public void test_getStartBasePos() {
		assertEquals(10, getStartBasePos(new BreakendSummary(0, FWD, 10, 20), 1));
		assertEquals(9, getStartBasePos(new BreakendSummary(0, FWD, 10, 20), 2));
	}
}

