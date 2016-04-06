package au.edu.wehi.idsv.picard;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class InMemoryReferenceSequenceFileTest extends TestHelper {
	@Test
	public void should_return_sequence() throws IOException {
		 try (InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(
				new String[] { "0", "1", },
				new byte[][] { B("CCCAATGGGCCC"),
							   B("TTTAATGGGAAA"), })) {
			assertEquals("CCCAATGGGCCC", S(ref.getSequence("0").getBases()));
			assertEquals("TTTAATGGGAAA", S(ref.getSequence("1").getBases()));
			assertEquals("C", S(ref.getSubsequenceAt("0", 1, 1).getBases()));
			assertEquals("C", S(ref.getSubsequenceAt("0", 2, 2).getBases()));
			assertEquals("C", S(ref.getSubsequenceAt("0", 3, 3).getBases()));
			assertEquals("A", S(ref.getSubsequenceAt("0", 4, 4).getBases()));
			assertEquals("CCCAAT", S(ref.getSubsequenceAt("0", 1, 6).getBases()));
		 }
	}
}
