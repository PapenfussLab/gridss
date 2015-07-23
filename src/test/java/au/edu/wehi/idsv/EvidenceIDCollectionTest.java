package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;


public class EvidenceIDCollectionTest extends TestHelper {
	@Test
	public void should_reverse_dp_sc() {
		EvidenceIDCollection col = new EvidenceIDCollection();
		EvidenceIDCollection rcol = new EvidenceIDCollection();
		DiscordantReadPair dp = (DiscordantReadPair)NRRP(SES(true), DP(0, 1, "1M", true, 0, 10, "1M", false));
		col.categorise(dp);
		rcol.categorise(dp.asRemote());
		rcol.categorise((DiscordantReadPair)NRRP(SES(false), DP(1, 1, "2M", true, 1, 10, "2M", false)));
		col.addRemote(rcol);
		assertArrayEquals(new int[] { 1, 1 }, col.getCount(VcfAttributes.BREAKPOINT_READPAIR_COUNT, 2));
	}
}
