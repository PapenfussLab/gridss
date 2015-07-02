package au.edu.wehi.idsv.picard;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.picard.BufferedReferenceSequenceFile;


public class BufferedReferenceSequenceFileTest extends TestHelper {
	@Test
	public void getSequenceShouldMatchUnderlying() throws IOException {
		BufferedReferenceSequenceFile b = new BufferedReferenceSequenceFile(SMALL_FA);
		assertEquals(S(b.getSequence("polyA").getBases()), S(SMALL_FA.getSequence("polyA").getBases()));
		assertEquals(S(b.getSequence("polyACGT").getBases()), S(SMALL_FA.getSequence("polyACGT").getBases()));
	}
	@Test
	public void getSubsequenceAtShouldMatchUnderlying() throws IOException {
		BufferedReferenceSequenceFile b = new BufferedReferenceSequenceFile(SMALL_FA);
		for (int i = 1; i < 100; i++) {
			for (int j = i; j < 100; j++) {
				assertEquals(S(b.getSubsequenceAt("polyACGT", i, j).getBases()), S(SMALL_FA.getSubsequenceAt("polyACGT", i, j).getBases()));
				assertEquals(b.getSubsequenceAt("polyACGT", i, j).getName(), SMALL_FA.getSubsequenceAt("polyACGT", i, j).getName());
				assertEquals(b.getSubsequenceAt("polyACGT", i, j).getContigIndex(), SMALL_FA.getSubsequenceAt("polyACGT", i, j).getContigIndex());
			}
		}
	}
}
