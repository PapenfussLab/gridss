package gridss.filter;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;


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
		assertFalse(new ClippedReadFilter(3, false).filterOut(Read(0, 1, "1H2S1M")));
		assertTrue(new ClippedReadFilter(3, false).filterOut(Read(0, 1, "2H1M")));
		assertTrue(new ClippedReadFilter(3, false).filterOut(Read(0, 1, "1H1M1S")));
		assertFalse(new ClippedReadFilter(3, false).filterOut(Read(0, 1, "2H1M1S")));
	}
	@Test
	public void should_filter_split_reads() {
		SAMRecord sc = Read(0, 1, "50M50S");
		SAMRecord sr = Read(0, 100, "50M50S");
		sr.setAttribute("SA", new ChimericAlignment(sc).toString());
		assertFalse(new ClippedReadFilter(3, true).filterOut(sc));
		assertFalse(new ClippedReadFilter(3, true).filterOut(sr));
		assertFalse(new ClippedReadFilter(3, false).filterOut(sc));
		assertTrue(new ClippedReadFilter(3, false).filterOut(sr));
	}
}
