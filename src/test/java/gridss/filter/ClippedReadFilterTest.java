package gridss.filter;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class ClippedReadFilterTest extends TestHelper {
	@Test
	public void should_filter_unclipped_read() {
		assertTrue(new ClippedReadFilter().filterOut(Read(0, 1, "1M")));
	}
	@Test
	public void should_not_filter_clipped_read() {
		assertFalse(new ClippedReadFilter().filterOut(Read(0, 1, "1S1M")));
		assertFalse(new ClippedReadFilter().filterOut(Read(0, 1, "1H1M")));
	}
	@Test
	public void should_filter_based_on_total_clipping_length() {
		assertFalse(new ClippedReadFilter(3).filterOut(Read(0, 1, "1H2S1M")));
		assertTrue(new ClippedReadFilter(3).filterOut(Read(0, 1, "2H1M")));
		assertTrue(new ClippedReadFilter(3).filterOut(Read(0, 1, "1H1M1S")));
		assertFalse(new ClippedReadFilter(3).filterOut(Read(0, 1, "2H1M1S")));
	}
}
