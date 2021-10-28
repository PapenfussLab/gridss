package au.edu.wehi.idsv.model;

import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class MapqModelTest {
	private static final MapqModel model = new MapqModel();
	@Test
	public void should_return_mapq() {
		assertEquals(5, model.scoreIndel(null, null, null, -1, 5), 0);
		assertEquals(5, model.scoreSoftClip(null,null, -1, 5), 0);
		assertEquals(5, model.scoreUnmappedMate(null,null, 5), 0);
	}
	@Test
	public void should_combine_mapq() {
		assertEquals(7.21246399, model.scoreSplitRead(null,null,  0, 10, 10), 0.000001);
		assertEquals(7.21246399, model.scoreReadPair(null,null,  0, 10, 10), 0.000001);
	}
}
