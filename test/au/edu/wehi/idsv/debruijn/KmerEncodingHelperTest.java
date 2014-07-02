package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.TestHelper;

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
		assertEquals(1,  KmerEncodingHelper.picardBaseToEncoded(1, new byte[] {'C'}));
		assertEquals(1,  KmerEncodingHelper.picardBaseToEncoded(1, new byte[] {'c'}));
	}
	@Test
	public void picardBaseToEncoded_should_initial_array_bases() {
		assertEquals(10,  KmerEncodingHelper.picardBaseToEncoded(2, B("AACC")));
	}
	@Test
	public void encoding_should_rountrip() {
		assertEquals("ACTG", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG")), 4)));
	}
	@Test
	public void nextState_should_clear_higher_bits() {
		assertEquals(15, KmerEncodingHelper.nextState(2, 15, (byte)'G'));
	}
	@Test
	public void nextState_should_handle_25mer() {
		assertEquals(750599937895082L, KmerEncodingHelper.nextState(25, 750599937895082L, (byte)'A'));
	}
	@Test
	public void nextState_should_handle_31mer() {
		assertEquals(KmerEncodingHelper.picardBaseToEncoded(31, B("TATATATATATATATATATATATATATATAG")), KmerEncodingHelper.nextState(31, KmerEncodingHelper.picardBaseToEncoded(31, B("ATATATATATATATATATATATATATATATA")), (byte)'G'));
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
	public void prev_state_should_advance_backward_3mer() {
		long state = KmerEncodingHelper.picardBaseToEncoded(3, B("AAC"));
		assertEquals("TAA", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 3)[0], 3)));
		assertEquals("CAA", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 3)[1], 3)));
		assertEquals("AAA", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 3)[2], 3)));
		assertEquals("GAA", S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.prevStates(state, 3)[3], 3)));
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
	@Test
	public void reverse_should_change_base_order() {
		String[] input =    new String[] {"AAAAAAAATTTTTTTTCCCCCCCCGGGGGGGG", "GTACGTACGTACGTACGTACGTACGTACGTAC", "CGCAGCTTATGTTGGGCGAGAACGCTAATTAC", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"};
		for (int i = 0; i < input.length; i++) {
			for (int k = 1; k < KmerEncodingHelper.MAX_KMER; k++) {
				assertEquals(new StringBuilder(input[i].substring(0, k)).reverse().toString(), S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.reverse(k, KmerEncodingHelper.picardBaseToEncoded(k, B(input[i].substring(0, k)))), k)));
				assertEquals(new StringBuilder(input[i].substring(k)).reverse().toString(), S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.reverse(KmerEncodingHelper.MAX_KMER-k, KmerEncodingHelper.picardBaseToEncoded(KmerEncodingHelper.MAX_KMER-k, B(input[i].substring(k)))), KmerEncodingHelper.MAX_KMER-k)));
			}
		}
	}
	@Test
	public void reverse_should_not_fail_on_sign_bit() {
		assertEquals(KmerEncodingHelper.picardBaseToEncoded(32, B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG")), KmerEncodingHelper.reverse(32, KmerEncodingHelper.picardBaseToEncoded(32, B("GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))));
		assertEquals(KmerEncodingHelper.picardBaseToEncoded(32, B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")), KmerEncodingHelper.reverse(32, KmerEncodingHelper.picardBaseToEncoded(32, B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))));
		assertEquals(KmerEncodingHelper.picardBaseToEncoded(32, B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT")), KmerEncodingHelper.reverse(32, KmerEncodingHelper.picardBaseToEncoded(32, B("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))));
		assertEquals(KmerEncodingHelper.picardBaseToEncoded(32, B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC")), KmerEncodingHelper.reverse(32, KmerEncodingHelper.picardBaseToEncoded(32, B("CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))));
	}
	@Test
	public void complement_should_change_bases() {
		String[] input =    new String[] {"CATTAATCGCAAGAGCGGGTTGTATTCGACGC", "GGGGGGGGCCCCCCCCTTTTTTTTAAAAAAAA", "CATGCATGCATGCATGCATGCATGCATGCATG", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"};
		String[] expected = new String[] {"GTAATTAGCGTTCTCGCCCAACATAAGCTGCG", "CCCCCCCCGGGGGGGGAAAAAAAATTTTTTTT", "GTACGTACGTACGTACGTACGTACGTACGTAC", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"};
		for (int i = 0; i < input.length; i++) {
			for (int k = 1; k < KmerEncodingHelper.MAX_KMER; k++) {
				assertEquals(expected[i].substring(0, k), S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.complement(k, KmerEncodingHelper.picardBaseToEncoded(k, B(input[i].substring(0, k)))), k)));
				assertEquals(expected[i].substring(k), S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.complement(KmerEncodingHelper.MAX_KMER-k, KmerEncodingHelper.picardBaseToEncoded(KmerEncodingHelper.MAX_KMER-k, B(input[i].substring(k)))),KmerEncodingHelper.MAX_KMER-k)));
			}
		}
	}
}
