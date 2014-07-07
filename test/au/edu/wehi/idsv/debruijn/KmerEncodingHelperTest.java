package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.apache.commons.lang3.StringUtils;
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
		assertEquals("ACTG", S(KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG")))));
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
		assertArrayEquals(B("CTGT"), KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.nextStates(4, state)[0]));
		assertArrayEquals(B("CTGC"), KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.nextStates(4, state)[1]));
		assertArrayEquals(B("CTGA"), KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.nextStates(4, state)[2]));
		assertArrayEquals(B("CTGG"), KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.nextStates(4, state)[3]));
	}
	@Test
	public void prev_state_should_advance_backward() {
		long state = KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG"));
		assertEquals("TACT", S(KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.prevStates(4, state)[0])));
		assertEquals("CACT", S(KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.prevStates(4, state)[1])));
		assertEquals("AACT", S(KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.prevStates(4, state)[2])));
		assertEquals("GACT", S(KmerEncodingHelper.encodedToPicardBases(4, KmerEncodingHelper.prevStates(4, state)[3])));
	}
	@Test
	public void prev_state_should_advance_backward_64bit() {
		long state = KmerEncodingHelper.picardBaseToEncoded(25, B("TACAATTCAGCATGAGCTAAACCCT"));
		assertEquals("TTACAATTCAGCATGAGCTAAACCC", S(KmerEncodingHelper.encodedToPicardBases(25, KmerEncodingHelper.prevStates(25, state)[0])));
		assertEquals("CTACAATTCAGCATGAGCTAAACCC", S(KmerEncodingHelper.encodedToPicardBases(25, KmerEncodingHelper.prevStates(25, state)[1])));
		assertEquals("ATACAATTCAGCATGAGCTAAACCC", S(KmerEncodingHelper.encodedToPicardBases(25, KmerEncodingHelper.prevStates(25, state)[2])));
		assertEquals("GTACAATTCAGCATGAGCTAAACCC", S(KmerEncodingHelper.encodedToPicardBases(25, KmerEncodingHelper.prevStates(25, state)[3])));
	}
	@Test
	public void prev_state_should_advance_backward_3mer() {
		long state = KmerEncodingHelper.picardBaseToEncoded(3, B("AAC"));
		assertEquals("TAA", S(KmerEncodingHelper.encodedToPicardBases(3, KmerEncodingHelper.prevStates(3, state)[0])));
		assertEquals("CAA", S(KmerEncodingHelper.encodedToPicardBases(3, KmerEncodingHelper.prevStates(3, state)[1])));
		assertEquals("AAA", S(KmerEncodingHelper.encodedToPicardBases(3, KmerEncodingHelper.prevStates(3, state)[2])));
		assertEquals("GAA", S(KmerEncodingHelper.encodedToPicardBases(3, KmerEncodingHelper.prevStates(3, state)[3])));
	}
	@Test
	public void adjacentStates_should_match_next_prevStates() {
		for (int k = 1; k <= 4; k++) {
			for (long state = 0; state < 2 << (2*k); state++) {
				assertEquals(8, KmerEncodingHelper.adjacentStates(k, state).length);
				assertEquals(KmerEncodingHelper.nextStates(k, state)[0], KmerEncodingHelper.adjacentStates(k, state)[0]);
				assertEquals(KmerEncodingHelper.nextStates(k, state)[1], KmerEncodingHelper.adjacentStates(k, state)[1]);
				assertEquals(KmerEncodingHelper.nextStates(k, state)[2], KmerEncodingHelper.adjacentStates(k, state)[2]);
				assertEquals(KmerEncodingHelper.nextStates(k, state)[3], KmerEncodingHelper.adjacentStates(k, state)[3]);
				assertEquals(KmerEncodingHelper.prevStates(k, state)[0], KmerEncodingHelper.adjacentStates(k, state)[4]);
				assertEquals(KmerEncodingHelper.prevStates(k, state)[1], KmerEncodingHelper.adjacentStates(k, state)[5]);
				assertEquals(KmerEncodingHelper.prevStates(k, state)[2], KmerEncodingHelper.adjacentStates(k, state)[6]);
				assertEquals(KmerEncodingHelper.prevStates(k, state)[3], KmerEncodingHelper.adjacentStates(k, state)[7]);
			}
		}
	}
	@Test
	public void firstBaseEncodedToPicardBase_should_return_first_base() {
		long state = KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG"));
		assertEquals((byte)'A', KmerEncodingHelper.firstBaseEncodedToPicardBase(4, state));
	}
	@Test
	public void lastBaseEncodedToPicardBase_should_return_first_base() {
		long state = KmerEncodingHelper.picardBaseToEncoded(4, B("ACTG"));
		assertEquals((byte)'G', KmerEncodingHelper.lastBaseEncodedToPicardBase(4, state));
	}
	@Test
	public void reverse_should_change_base_order() {
		String[] input =    new String[] {"AAAAAAAATTTTTTTTCCCCCCCCGGGGGGGG", "GTACGTACGTACGTACGTACGTACGTACGTAC", "CGCAGCTTATGTTGGGCGAGAACGCTAATTAC", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"};
		for (int i = 0; i < input.length; i++) {
			for (int k = 1; k < KmerEncodingHelper.MAX_K; k++) {
				assertEquals(new StringBuilder(input[i].substring(0, k)).reverse().toString(), S(KmerEncodingHelper.encodedToPicardBases(k, KmerEncodingHelper.reverse(k, KmerEncodingHelper.picardBaseToEncoded(k, B(input[i].substring(0, k)))))));
				assertEquals(new StringBuilder(input[i].substring(k)).reverse().toString(), S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.MAX_K-k, KmerEncodingHelper.reverse(KmerEncodingHelper.MAX_K-k, KmerEncodingHelper.picardBaseToEncoded(KmerEncodingHelper.MAX_K-k, B(input[i].substring(k)))))));
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
			for (int k = 1; k < KmerEncodingHelper.MAX_K; k++) {
				assertEquals(expected[i].substring(0, k), S(KmerEncodingHelper.encodedToPicardBases(k, KmerEncodingHelper.complement(k, KmerEncodingHelper.picardBaseToEncoded(k, B(input[i].substring(0, k)))))));
				assertEquals(expected[i].substring(k), S(KmerEncodingHelper.encodedToPicardBases(KmerEncodingHelper.MAX_K-k,KmerEncodingHelper.complement(KmerEncodingHelper.MAX_K-k, KmerEncodingHelper.picardBaseToEncoded(KmerEncodingHelper.MAX_K-k, B(input[i].substring(k)))))));
			}
		}
	}
	public void assertBasesMatching(int expectedMismatches, String a, String b) {
		assertEquals(expectedMismatches, KmerEncodingHelper.basesDifference(a.length(),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(a)),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(b))));
		assertEquals(a.length() - expectedMismatches, KmerEncodingHelper.basesMatching(a.length(),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(a)),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(b))));
		assertEquals(expectedMismatches, KmerEncodingHelper.basesDifference(a.length(),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(b)),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(a))));
		assertEquals(a.length() - expectedMismatches, KmerEncodingHelper.basesMatching(a.length(),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(b)),
				KmerEncodingHelper.picardBaseToEncoded(a.length(), B(a))));
	}
	@Test
	public void basesDifference_basesMatching_should_work_at_base_level() {
		assertBasesMatching(0, "A", "A");
		assertBasesMatching(1, "G", "A");
		assertBasesMatching(0, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAATAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		
		assertBasesMatching(0, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAATAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAGAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		assertBasesMatching(1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		
		assertBasesMatching(31, "AAAAAATAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(30, "TAAAAATAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(30, "AAAAAATAAAAAAAAAATAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(30, "AAAAAATAAAAAAAAAAAAAAAAAAAAAAAAT", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(31, "GAAAAATAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(31, "AAAAAATAAAAAAAAAAGAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(31, "AAAAAATAAAAAAAAAAAAAAAAAAAAAAAAG", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(31, "CAAAAATAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(31, "ACAAAATAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertBasesMatching(31, "AAAAAATAAAAAAAAAAAAAAAAAAAAAAAAC", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		
		assertBasesMatching(12, "GTACGGTACAGTACTGTACC", "GGGGGAAAAATTTTTCCCCC");
	}
}
