package scambler;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class StringGraphTransitiveCompressorTest extends StringGraphTestHelper {
	@Test
	public void should_return_all_nodes() {
		for (int i = 2; i < 10; i++) {
			StringGraphTransitiveCompressor trit = new StringGraphTransitiveCompressor(
					new StreamingStringGraphIterator(10, 200, 30, 
						overlapping(i, 30, 10).iterator(),
						getContext().getLinear()),
					0, false);
			ArrayList<SgNode> result = Lists.newArrayList(trit);
			assertEquals(2 * i, result.size());
		}
	}
	@Test
	public void should_transitive_reduce_simple_three_read_overlap() {
		// AAA
		//  BBB
		//   CCC
		List<SAMRecord> reads = overlapping(3, 9, 3);
		ArrayList<SgNode> rawGraph = Lists.newArrayList(new StreamingStringGraphIterator(3, 200, 3, reads.iterator(), getContext().getLinear()));
		StringGraphTransitiveCompressor sit = new StringGraphTransitiveCompressor(rawGraph.iterator(), 0, false);
		ArrayList<SgNode> graph = Lists.newArrayList(sit);
		assertEquals(6, graph.size());
		graph.sort(SgNode.ByInferredPosition);
		SgNode aStart = graph.get(0);
		SgNode bStart = graph.get(1);
		SgNode cStart = graph.get(2);
		SgNode aEnd = graph.get(3);
		SgNode bEnd = graph.get(4);
		SgNode cEnd = graph.get(5);
		assertEquals(bStart, aStart.out.get(0).to);
		//assertEquals(cStart, aStart.out.get(1).to);
		assertEquals(cStart, bStart.out.get(0).to);
		//assertEquals(aEnd, bStart.out.get(1).to);
		assertEquals(aEnd, cStart.out.get(0).to);
		//assertEquals(bEnd, cStart.out.get(1).to);
		assertEquals(bEnd, aEnd.out.get(0).to);
		//assertEquals(cEnd, aEnd.out.get(1).to);
		assertEquals(cEnd, bEnd.out.get(0).to);
		assertEquals(5, graph.stream().mapToInt(n -> n.out.size()).sum());
	}
	@Test
	public void should_pass_through_simple_two_read_overlap() {
		List<SAMRecord> reads = overlapping(2, 9, 3);
		ArrayList<SgNode> rawGraph = Lists.newArrayList(new StreamingStringGraphIterator(3, 200, 3, reads.iterator(), getContext().getLinear()));
		// sanity check graph construction
		assertEquals(4, rawGraph.size());
		assertEquals(3, rawGraph.stream().mapToInt(n -> n.out.size()).sum());
		StringGraphTransitiveCompressor tcit = new StringGraphTransitiveCompressor(rawGraph.iterator(), 0, false);
		ArrayList<SgNode> graph = Lists.newArrayList(tcit);
		assertEquals(4, graph.size());
		assertEquals(3, graph.stream().mapToInt(n -> n.out.size()).sum());
	}
	@Test
	@Ignore("NYI")
	public void should_reduce_perfect_graph_to_single_path() {
		for (int i = 3; i < 100; i++) {
			List<SAMRecord> reads = overlapping(i, 60, 20);
			ArrayList<SgNode> rawGraph = Lists.newArrayList(new StreamingStringGraphIterator(20, 60, 60, reads.iterator(), getContext().getLinear()));
			StringGraphTransitiveCompressor tcit = new StringGraphTransitiveCompressor(rawGraph.iterator(), 0, false);
			ArrayList<SgNode> graph = Lists.newArrayList(tcit);
			assertEquals(graph.size() - 1, graph.stream().mapToInt(n -> n.out.size()).sum());
		}
	}
	@Test
	@Ignore("NYI")
	public void should_reduce_perfect_graph_to_single_node() {
		for (int i = 3; i < 100; i++) {
			List<SAMRecord> reads = overlapping(i, 60, 20);
			ArrayList<SgNode> rawGraph = Lists.newArrayList(new StreamingStringGraphIterator(20, 60, 60, reads.iterator(), getContext().getLinear()));
			StringGraphTransitiveCompressor tcit = new StringGraphTransitiveCompressor(rawGraph.iterator(), 0, true);
			ArrayList<SgNode> graph = Lists.newArrayList(tcit);
			assertEquals(0, graph.size());
		}
	}
}
