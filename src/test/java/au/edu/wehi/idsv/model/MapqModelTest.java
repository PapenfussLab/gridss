package au.edu.wehi.idsv.model;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class MapqModelTest {
	private static final MapqModel model = new MapqModel();
	@Test
	public void should_return_mapq() {
		assertEquals(5, model.scoreIndel(null, null, -1, 5), 0);
		assertEquals(5, model.scoreSoftClip(null, -1, 5), 0);
		assertEquals(5, model.scoreUnmappedMate(null, 5), 0);
	}
	@Test
	public void should_combine_mapq() {
		assertEquals(7.21246399, model.scoreSplitRead(null, 0, 10, 10), 0.000001);
		assertEquals(7.21246399, model.scoreReadPair(null, 0, 10, 10), 0.000001);
	}
}
