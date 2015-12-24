package au.edu.wehi.idsv.model;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class MapqModelTest {
	@Test
	public void should_return_mapq() {
		assertEquals(5, new MapqModel().scoreIndel(null, null, -1, 5), 0.0001);
	}
	@Test
	public void should_combine_mapq() {
		assertEquals(7.21246399, new MapqModel().scoreSplitRead(null, 0, 10, 10), 0.0001);
	}
}
