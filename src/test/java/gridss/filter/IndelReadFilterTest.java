package gridss.filter;

import au.edu.wehi.idsv.TestHelper;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;


public class IndelReadFilterTest extends TestHelper {
	@Test
	public void should_filter_unclipped_read() {
		assertTrue(new IndelReadFilter().filterOut(Read(0, 1, "1M")));
	}
	@Test
	public void should_not_filter_indel_reads() {
		assertFalse(new IndelReadFilter().filterOut(Read(0, 1, "1M1I1M")));
		assertFalse(new IndelReadFilter().filterOut(Read(0, 1, "1M1D1M")));
	}
	@Test
	public void should_not_filter_skip_reads() {
		assertFalse(new IndelReadFilter().filterOut(Read(0, 1, "1M1N1M")));
	}
	@Test
	public void should_filter_based_on_max_indel_length() {
		assertTrue(new IndelReadFilter(2).filterOut(Read(0, 1, "1M1I1D1M")));
		assertFalse(new IndelReadFilter(2).filterOut(Read(0, 1, "1M2I1M")));
	}
}