package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;


public class PackedSequenceTest extends TestHelper {
	@Test
	public void should_round_trip_sequence() {
		for (String s : KmerEncodingHelperTest.TWOMERS) {
			PackedSequence seq = new PackedSequence(B(s), false, false);
			assertEquals(s, S(seq.getBytes(0, s.length())));
			assertEquals(s.charAt(0), (char)seq.get(0));
			assertEquals(s.charAt(1), (char)seq.get(1));
		}
		for (String s : ImmutableList.of(
				"A", "C", "G", "T",
				S(POLY_A),
				S(POLY_ACGT),
				S(RANDOM))) {
			PackedSequence seq = new PackedSequence(B(s), false, false);
			assertEquals(s, S(seq.getBytes(0, s.length())));
			assertEquals(s.substring(1), S(seq.getBytes(1, s.length() - 1)));
		}
	}
	@Test
	public void should_allow_1_to_32_base_kmers() {
		String seq = "CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGATTTTGTTTACAGCCTGTCTTATATCCTGAATAACGCACCGCCTATTCGAACGGGCGAATCTACCTAGGTCGCTCAGAACCGGCACCCTTAACCATCCATAT";
		PackedSequence list = new PackedSequence(B(seq), false, false);
		for (int k = 1; k <= 32; k++) {
			for (int i = 0; i < seq.length()-(k-1); i++) {
				long kmer = list.getKmer(i, k);
				assertEquals(KmerEncodingHelper.picardBaseToEncoded(k, B(seq.substring(i, i+k))), kmer);
			}
		}
	}
	@Test
	public void should_reverse() {
		PackedSequence seq = new PackedSequence(B("AACGT"), true, false);
		assertEquals("TGCAA", S(seq.getBytes(0, 5)));
	}
	@Test
	public void should_comp() {
		PackedSequence seq = new PackedSequence(B("AACGT"), false, true);
		assertEquals("TTGCA", S(seq.getBytes(0, 5)));
	}
	@Test
	public void should_allow_greater_than_word_bases() {
		PackedSequence seq = new PackedSequence(B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC"), false, false);
		assertEquals("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC", S(seq.getBytes(0, 33)));
	}
	@Test
	public void overlapMatches_should_return_matching_base_count() {
		PackedSequence seq1 = new PackedSequence(B("ACGT"), false, false);
		assertEquals(4, PackedSequence.overlapMatches(seq1, seq1, 0));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, 1));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, 2));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, 3));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, 4));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, -1));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, -2));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, -3));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq1, -4));
	}
	@Test
	public void overlapMatches_should_allow_difference_sequence_lengths() {
		PackedSequence seq1 = new PackedSequence(B("ACGT"), false, false);
		PackedSequence seq2 = new PackedSequence(B("ACTTT"), false, false);
		// ACGT
		// ACTTT
		assertEquals(3, PackedSequence.overlapMatches(seq1, seq2, 0));
		// ACGT
		// -ACTTT
		assertEquals(1, PackedSequence.overlapMatches(seq1, seq2, 1));
		// ACGT
		// --ACTTT
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, 2));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, 3));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, 4));
		// -ACGT
		// ACTTT
		assertEquals(1, PackedSequence.overlapMatches(seq1, seq2, -1));
		// --ACGT
		// ACTTT
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, -2));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, -3));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, -4));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, -5));
		assertEquals(0, PackedSequence.overlapMatches(seq1, seq2, -6));
	}
	@Test
	public void overlapMatches_should_span_words() {
		PackedSequence seq1 = new PackedSequence(B(
				"CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGATTTTGTTTACAGCCTGTCTTATATCCTGAATAACGCACCGCCTATTCGAACGGGCGAATCTACCTAGGTCGCTCAGAACCGGCACCCTTAACCATCCATAT"), false, false);
		PackedSequence seq2 = new PackedSequence(B(
				"CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAAAAAGTTTACAGCCTGTCTTATATCCTGAATAACGCACCGCCTATTCGAACGGGCGAATCTACCTAGGTCGCTCAGAACCGGCACCCTTAACCATCCATAT"), false, false);
		assertEquals(seq1.length(), PackedSequence.overlapMatches(seq1, seq1, 0));
		assertEquals(seq1.length() - 4, PackedSequence.overlapMatches(seq1, seq2, 0));
	}
	@Test
	public void overlapMatches_should_match_offset_reads() {
		PackedSequence seq1 = new PackedSequence(B("TTATTTGAGTTAGAGAGCTTCGGGGATACTTATTCAGGTACGGTCTGTCTTCGGGAGTGTGACAGTACATGAATTAGCCTTGGTAGCAAGACAATTTTGT"), false, false);
		PackedSequence seq2 = new PackedSequence(B( "TATTTGAGTTAGAGAGCTTCGGGGATACTTATTCAGGTACGGTCTGTCTTCGGGAGTGTGACAGTACATGAATTAGCCTTGGTAGCAAGACAATTTTGTG"), false, false);
		assertEquals(99, PackedSequence.overlapMatches(seq1, seq2, 1));
	}
	private void test_overlapLength(String s1, String s2, int offset, int expectedOverlap) {
		PackedSequence seq1 = new PackedSequence(B(s1), false, false);
		PackedSequence seq2 = new PackedSequence(B(s2), false, false);
		assertEquals(expectedOverlap, PackedSequence.overlapLength(seq1, seq2, offset));
	}
	@Test
	public void overlapLength_should_return_common_base_count() {
		test_overlapLength("AAAT", "AA", -3, 0);
		test_overlapLength("AAAT", "AA", -2, 0);
		test_overlapLength("AAAT", "AA", -1, 1);
		test_overlapLength("AAAT", "AA", 0, 2);
		test_overlapLength("AAAT", "AA", 1, 2);
		test_overlapLength("AAAT", "AA", 2, 2);
		test_overlapLength("AAAT", "AA", 3, 1);
		test_overlapLength("AAAT", "AA", 4, 0);
		test_overlapLength("AAAT", "AA", 5, 0);
	}
	@Test
	public void subsequence_constructor_should_return_subsequence() {
		for (int len = 0; len < 100; len++) {
			PackedSequence seq = new PackedSequence(B(S(RANDOM).substring(0, len)), false, false);
			for (int i = 0; i < len; i++) {
				for (int j = 0; i + j <= len; j++) {
					PackedSequence subseq = new PackedSequence(seq, i, j);
					assertEquals(seq.toString().substring(i, i + j), subseq.toString());
				}
			}
		}
	}
	@Test
	public void should_concat_sequences() {
		for (int i =  0; i < 200; i++) {
			for (int j =  0; j < 200; j++) {
				String seq1 = S(RANDOM).substring(i, i + j);
				String seq2 = S(RANDOM).substring(j, i + j);
				PackedSequence ps1 = new PackedSequence(B(seq1), false, false);
				PackedSequence ps2 = new PackedSequence(B(seq2), false, false);
				PackedSequence concat = new PackedSequence(ps1, ps2);
				assertEquals(seq1 + seq2, S(concat.getBytes(0, concat.length())));
			}
		}
	}
}
