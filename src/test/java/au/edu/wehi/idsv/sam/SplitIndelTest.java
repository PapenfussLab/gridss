package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class SplitIndelTest extends TestHelper {
	@Test
	public void shouldSplitAtDeletion() {
		List<SplitIndel> list = SplitIndel.getIndelsAsSplitReads(Read(0, 1, "1M1D1M"));
		assertEquals("1M1S", list.get(0).leftAnchored.getCigarString());
		assertEquals("1S1M", list.get(0).rightAnchored.getCigarString());
		assertEquals(1, list.get(0).leftAnchored.getAlignmentStart());
		assertEquals(3, list.get(0).leftRealigned.getAlignmentStart());
		assertEquals(1, list.get(0).rightRealigned.getAlignmentStart());
		assertEquals(3, list.get(0).rightAnchored.getAlignmentStart());
	}
	@Test
	public void shouldSplitAtInsertion() {
		List<SplitIndel> list = SplitIndel.getIndelsAsSplitReads(Read(0, 1, "1M1I1M"));
		assertEquals("1M2S", list.get(0).leftAnchored.getCigarString());
		assertEquals("2S1M", list.get(0).rightAnchored.getCigarString());
		assertEquals(1, list.get(0).leftAnchored.getAlignmentStart());
		assertEquals(2, list.get(0).leftRealigned.getAlignmentStart());
		assertEquals(1, list.get(0).rightRealigned.getAlignmentStart());
		assertEquals(2, list.get(0).rightAnchored.getAlignmentStart());
	}
	@Test
	public void shouldOffsetRightAnchorAlignmentStart() {
		List<SplitIndel> list = SplitIndel.getIndelsAsSplitReads(Read(0, 100, "1M2P2N2P1M"));
		assertEquals(99, list.get(0).rightAnchored.getAlignmentStart());
		assertEquals(99, list.get(0).leftRealigned.getAlignmentStart());
	}
	@Test
	public void should_decode_negative_deletion() {
		List<SplitIndel> list = SplitIndel.getIndelsAsSplitReads(Read(0, 1000, "136M8952P8952N8952P136M"));
		assertEquals(1, list.size());
	}
	@Test
	public void should_merge_adjacent_indels() {
		List<SplitIndel> list = SplitIndel.getIndelsAsSplitReads(Read(0, 1000, "136M88I8952P8952N8952P10I136M"));
		assertEquals(1, list.size());
	}
	@Test
	public void should_split_small_insertion() {
		List<SplitIndel> list = SplitIndel.getIndelsAsSplitReads(Read(0, 1, "2M2I2M"));
		assertEquals(1, list.size());
		assertEquals("2M4S", list.get(0).leftAnchored.getCigarString());
		assertEquals("2S2M", list.get(0).leftRealigned.getCigarString());
		assertEquals("4S2M", list.get(0).rightAnchored.getCigarString());
		assertEquals("2M2S", list.get(0).rightRealigned.getCigarString());
	}
}
