package au.edu.wehi.idsv;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;


public class DiscordantReadPairTest extends TestHelper {
	public DiscordantReadPair newPair(SAMRecord[] pair, int maxfragmentSize) {
		return (DiscordantReadPair)NonReferenceReadPair.create(pair[0], pair[1], SES(maxfragmentSize));
	}
	@Test
	public void getRemoteMapq_should_be_anchored_mapq() {
		SAMRecord[] pair = DP(1, 1, "100M", true, 2, 5, "100M", true);
		pair[0].setMappingQuality(5);
		pair[1].setMappingQuality(10);
		assertEquals(10, newPair(pair, 300).getRemoteMapq());
	}
	@Test
	public void getRemoteBaseLength_should_be_read_length() {
		assertEquals(100, newPair(DP(1, 1, "50M", true, 2, 5, "100M", true), 300).getRemoteBaseLength());
	}
	@Test
	public void getRemoteBaseCount_should_be_read_length() {
		assertEquals(100, newPair(DP(1, 1, "100M", true, 2, 5, "100M", true), 300).getRemoteBaseCount());
	}
	@Test
	public void getRemoteMaxBaseQual_local_mapped_quals() {
		SAMRecord[] pair = DP(1, 1, "3M1S", true, 2, 5, "3M1S", true);
		withQual(new byte[] { 1, 2, 3, 4}, pair[0]);
		withQual(new byte[] { 4, 6, 5, 7}, pair[1]);
		assertEquals(6, newPair(pair, 300).getRemoteMaxBaseQual());
	}
	@Test
	public void getRemoteTotalBaseQual_local_mapped_quals() {
		SAMRecord[] pair = DP(1, 1, "3M1S", true, 2, 5, "1S3M", true);
		withQual(new byte[] { 1, 2, 3, 4}, pair[0]);
		withQual(new byte[] { 4, 5, 6, 7}, pair[1]);
		assertEquals(5+6+7, newPair(pair, 300).getRemoteTotalBaseQual());
	}
	@Test
	public void fragmentSequencesOverlap() {
		assertTrue(newPair(DP(1, 1, "5M", true, 1, 5, "5M", false), 300).fragmentSequencesOverlap());
		assertTrue(newPair(DP(1, 1, "5M", true, 1, 5, "5M", true), 300).fragmentSequencesOverlap());
		assertFalse(newPair(DP(1, 1, "5M", true, 0, 5, "5M", false), 300).fragmentSequencesOverlap());
		assertFalse(newPair(DP(1, 1, "5M", true, 1, 6, "5M", true), 300).fragmentSequencesOverlap());
	}
}
