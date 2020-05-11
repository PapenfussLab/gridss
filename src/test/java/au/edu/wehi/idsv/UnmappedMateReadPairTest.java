package au.edu.wehi.idsv;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class UnmappedMateReadPairTest extends TestHelper {
	@Test
	public void should_clip_breakend_to_contig_bounds() {
		UnmappedMateReadPair rp = (UnmappedMateReadPair) NRRP(OEA(0, 1, "1M", false));
		assertEquals(new BreakendSummary(0, BWD, 1), rp.getBreakendSummary());
	}
}
