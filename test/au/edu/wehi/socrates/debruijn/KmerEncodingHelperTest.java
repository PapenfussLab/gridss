package au.edu.wehi.socrates.debruijn;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.socrates.TestHelper;

public class KmerEncodingHelperTest extends TestHelper {
	@Test
	public void picardBaseToEncoded_should_use_blat_2bit_base_encoding() {
		assertEquals(0,  KmerEncodingHelper.picardBaseToEncoded(1, new byte[] {'T'}));
		assertEquals(1,  KmerEncodingHelper.picardBaseToEncoded(1, new byte[] {'C'}));
		assertEquals(2,  KmerEncodingHelper.picardBaseToEncoded(1, new byte[] {'A'}));
		assertEquals(3,  KmerEncodingHelper.picardBaseToEncoded(1, new byte[] {'G'}));
	}
	@Test
	public void picardBaseToEncoded_should_allow_both_upper_and_lower_case() {
		assertEquals(0,  KmerEncodingHelper.picardBaseToEncoded(1, new byte[] {'T'}));
	}
	@Test
	public void encoding_should_rountrip() {
		assertEquals("ACTG", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG")), 4)));
	}
	@Test
	public void next_state_should_advance_forward() {
		long state = KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG"));
		assertArrayEquals(B("CTGT"), KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.nextStates(state, 4)[0], 4));
		assertArrayEquals(B("CTGC"), KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.nextStates(state, 4)[1], 4));
		assertArrayEquals(B("CTGA"), KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.nextStates(state, 4)[2], 4));
		assertArrayEquals(B("CTGG"), KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.nextStates(state, 4)[3], 4));
	}
	@Test
	public void prev_state_should_advance_backward() {
		long state = KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG"));
		assertEquals("TACT", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 4)[0], 4)));
		assertEquals("CACT", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 4)[1], 4)));
		assertEquals("AACT", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 4)[2], 4)));
		assertEquals("GACT", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 4)[3], 4)));
	}
	@Test
	public void firstBaseEncodedToPicardBase_should_return_first_base() {
		long state = KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG"));
		assertEquals((byte)'A', KmerEncodingHelper.firstBaseEncodedToPicardBase(state, 4));
	}
	@Test
	public void lastBaseEncodedToPicardBase_should_return_first_base() {
		long state = KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG"));
		assertEquals((byte)'G', KmerEncodingHelper.lastBaseEncodedToPicardBase(state, 4));
	}
}
