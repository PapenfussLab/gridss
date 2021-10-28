package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class ReadCountModelTest {
	private static final ReadCountModel model = new ReadCountModel();
	@Test
	public void should_return_read_count() {
		assertEquals(1, model.scoreIndel(null, null, null, -1, 5), 0);
		assertEquals(1, model.scoreSoftClip(null,null, -1, 5), 0);
		assertEquals(1, model.scoreUnmappedMate(null,null, 5), 0);
		assertEquals(1, model.scoreSplitRead(null,null,  0, 10, 10), 0.000001);
		assertEquals(1, model.scoreReadPair(null,null,  0, 10, 10), 0.000001);
		assertEquals(5, model.scoreAssembly(null, 1, 2, 4, 8, 16, 32), 0.000001);
		assertEquals(5, model.scoreBreakendAssembly(null, 1, 2, 4, 8, 16), 0.000001);
	}
}
