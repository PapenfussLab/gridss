package scambler;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class StreamingStringGraphIteratorTest extends StringGraphTestHelper {
	@Test
	public void should_generate_two_nodes_per_read() {
		List<SAMRecord> reads = overlapping(RANDOM, 50, 100, 1, false);
		StreamingStringGraphIterator it = new StreamingStringGraphIterator(16, 200, 50, reads.iterator(), getContext().getLinear());
		ArrayList<SgNode> result = Lists.newArrayList(it);
		assertEquals(result.size(), 2 * 100);
	}
	@Test
	public void nodes_should_be_ordered_by_inferred_position() {
		List<SAMRecord> reads = overlapping(RANDOM, 50, 100, 1, false);
		StreamingStringGraphIterator it = new StreamingStringGraphIterator(16, 200, 50, reads.iterator(), getContext().getLinear());
		ArrayList<SgNode> result = Lists.newArrayList(it);
		assertTrue(SgNode.ByInferredPosition.isOrdered(result));
	}
	@Test
	public void edges_should_be_created() {
		List<SAMRecord> reads = overlapping(RANDOM, 50, 100, 1, true);
		StreamingStringGraphIterator it = new StreamingStringGraphIterator(48, 1000000000, 50, reads.iterator(), getContext().getLinear());
		ArrayList<SgNode> result = Lists.newArrayList(it);
		// 100 nodes * (before, common, after minus) * (two successor nodes each with 48 base overlaps)
		assertEquals(100 * 3 * 2, result.stream().mapToInt(node -> node.out.size()).sum());
	}
	@Test
	public void should_generate_edges_for_overlaps() {
		// AAA
		//  BBB
		//   CCC
		List<SAMRecord> reads = overlapping(3, 9, 3);
		StreamingStringGraphIterator sit = new StreamingStringGraphIterator(3, 200, 3, reads.iterator(), getContext().getLinear());
		ArrayList<SgNode> rawGraph = Lists.newArrayList(sit);
		assertEquals(6, rawGraph.size());
		rawGraph.sort(SgNode.ByInferredPosition);
		rawGraph.stream().forEach(n -> n.out.sort(SgEdge.BySequenceLength));
		SgNode aStart = rawGraph.get(0);
		SgNode bStart = rawGraph.get(1);
		SgNode cStart = rawGraph.get(2);
		SgNode aEnd = rawGraph.get(3);
		SgNode bEnd = rawGraph.get(4);
		SgNode cEnd = rawGraph.get(5);
		assertEquals(bStart, aStart.out.get(0).to);
		assertEquals(cStart, aStart.out.get(1).to);
		assertEquals(cStart, bStart.out.get(0).to);
		assertEquals(aEnd, bStart.out.get(1).to);
		assertEquals(aEnd, cStart.out.get(0).to);
		assertEquals(bEnd, cStart.out.get(1).to);
		assertEquals(bEnd, aEnd.out.get(0).to);
		assertEquals(cEnd, aEnd.out.get(1).to);
		assertEquals(cEnd, bEnd.out.get(0).to);
		assertEquals(9, rawGraph.stream().mapToInt(n -> n.out.size()).sum());
	}
	@Test
	public void perfect_overlaps_should_follow_expected_count() {
		for (int i = 3; i < 128; i++) {
			List<SAMRecord> reads = overlapping(i, 60, 20);
			StreamingStringGraphIterator sit = new StreamingStringGraphIterator(20, 100000, 60, reads.iterator(), getContext().getLinear());
			ArrayList<SgNode> rawGraph = Lists.newArrayList(sit);
			assertEquals(2 * i, rawGraph.size());
			// each read except the last two has two overlapping successors
			int overlaps = 2 * i - 3;
			// each overlap results in 3 edges
			assertEquals(3 * overlaps, rawGraph.stream().mapToInt(n -> n.out.size()).sum());
		}
	}
	@Test
	public void windowing_should_not_drop_edges() {
		List<SAMRecord> reads = overlapping(100, 60, 20);
		StreamingStringGraphIterator sit = new StreamingStringGraphIterator(20, 100000, 60, reads.iterator(), getContext().getLinear());
		ArrayList<SgNode> rawGraph = Lists.newArrayList(sit);
		
		List<SAMRecord> reads2 = overlapping(100, 60, 20);
		StreamingStringGraphIterator sit2 = new StreamingStringGraphIterator(20, 60, 60, reads2.iterator(), getContext().getLinear());
		ArrayList<SgNode> rawGraph2 = Lists.newArrayList(sit2);
		
		assertEquals(rawGraph.size(), rawGraph2.size());
		assertEquals(rawGraph.stream().mapToInt(n -> n.out.size()).sum(), rawGraph2.stream().mapToInt(n -> n.out.size()).sum());
	}
	@Test
	public void should_remove_non_overlapping_reads() {
		List<SAMRecord> reads = overlapping(10, 50, 50);
		StreamingStringGraphIterator sit = new StreamingStringGraphIterator(50, 100000, 50, reads.iterator(), getContext().getLinear());
		ArrayList<SgNode> rawGraph = Lists.newArrayList(sit);
		assertEquals(0, rawGraph.size());
	}
}
