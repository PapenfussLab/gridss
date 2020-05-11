package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class DiscordantReadPairTest extends TestHelper {
	public DiscordantReadPair newPair(SAMRecord[] pair, int maxfragmentSize) {
		return (DiscordantReadPair)NonReferenceReadPair.create(pair[0], pair[1], SES(maxfragmentSize));
	}
	@Test
	public void getRemoteMapq_should_be_anchored_mapq() {
		SAMRecord[] pair = DP(1, 1, "100M", true, 2, 5, "100M", true);
		pair[0].setMappingQuality(15);
		pair[1].setMappingQuality(20);
		assertEquals(20, newPair(pair, 300).getRemoteMapq());
	}
	@Test
	public void should_clip_breakpoint_to_contig_bounds() {
		SAMRecord[] pair = DP(0, 1, "1M", false, 0, 10000, "1M", true);
		DiscordantReadPair rp = newPair(pair, 300);
		assertEquals(new BreakpointSummary(0, BWD, 1, 0, FWD, 10000), rp.getBreakendSummary());
	}
	@Test
	public void remoteEvidenceID_should_evidenceID_of_other_side() {
		SAMRecord[] pair = DP(0, 1, "1M", false, 0, 10000, "1M", true);
		DiscordantReadPair rpf = (DiscordantReadPair)NonReferenceReadPair.create(pair[0], pair[1], SES(100));
		DiscordantReadPair rpb = (DiscordantReadPair)NonReferenceReadPair.create(pair[1], pair[0], SES(100));
		assertEquals(rpb.getEvidenceID(), rpf.getRemoteEvidenceID());
		assertEquals(rpf.getEvidenceID(), rpb.getRemoteEvidenceID());
	}
}
