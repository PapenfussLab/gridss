package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class SequenceUtilTest extends TestHelper {
	public static final double DELTA = 1e-9; // differs to excel calculation at 1e-10 level 
	@Test
	public void shannonEntropy() {
		assertEquals(0, SequenceUtil.shannonEntropy(B(""), 0, 0), DELTA);
		assertEquals(1.846439345, SequenceUtil.shannonEntropy(B("AGGCCCTTTT"), 0, 10), DELTA);
		assertEquals(1.846439345, SequenceUtil.shannonEntropy(B("AGGCCCTTTTAAAAA"), 0, 10), DELTA);
		assertEquals(1.846439345, SequenceUtil.shannonEntropy(B("TAGGCCCTTTTAAAAA"), 1, 10), DELTA);
		assertEquals(2, SequenceUtil.shannonEntropy(B("ACGT"), 0, 4), DELTA);
		assertEquals(1.921928095, SequenceUtil.shannonEntropy(B("ACTGT"), 0, 5), DELTA);
		assertEquals(1.921928095, SequenceUtil.shannonEntropy(B("TACGT"), 0, 5), DELTA);
	}
}
