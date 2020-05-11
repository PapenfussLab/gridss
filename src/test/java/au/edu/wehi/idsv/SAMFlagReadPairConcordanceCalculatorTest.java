package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

public class SAMFlagReadPairConcordanceCalculatorTest extends TestHelper {
	@Test
	public void should_require_same_reference_sequence() {
		SAMRecord[] rp = RP(0, 10, 100);
		rp[0].setProperPairFlag(true);
		Assert.assertTrue(new SAMFlagReadPairConcordanceCalculator(null).isConcordant(rp[0]));
	}
}
