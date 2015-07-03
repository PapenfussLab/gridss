package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.*;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;


public class PackedSequenceTest extends TestHelper {
	@Test
	public void should_round_trip_sequence() {
		for (String s : KmerEncodingHelperTest.TWOMERS) {
			PackedSequence seq = new PackedSequence(B(s), false, false);
			assertEquals(s, S(seq.getBytes(0, s.length())));
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
}
