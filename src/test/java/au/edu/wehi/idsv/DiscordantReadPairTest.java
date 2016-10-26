package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;


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
}
