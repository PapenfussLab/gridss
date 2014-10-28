package au.edu.wehi.idsv;


import static org.junit.Assert.*;

import org.junit.Test;


public class UnmappedMateReadPairTest extends TestHelper {
	@Test
	public void fragmentSequencesOverlap() {
		assertFalse(new UnmappedMateReadPair(OEA(0, 1, "5M", true)[0], OEA(0, 1, "5M", true)[1], SES(300)).fragmentSequencesOverlap());
	}
}
