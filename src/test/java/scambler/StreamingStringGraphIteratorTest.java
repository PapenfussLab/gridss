package scambler;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMRecord;

public class StreamingStringGraphIteratorTest extends TestHelper {
	private List<SAMRecord> overlapping(byte[] seq, int readLength, int seqlength, int stride, boolean isCircular) {
		String seq1 = S(seq).substring(0, seqlength);
		String seq2 = S(seq).substring(seqlength, 2 * seqlength);
		if (isCircular) {
			seq2 = seq1;
		}
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < seq1.length(); i += stride) {
			String perfect_read = (seq1 + seq2).substring(i, i + readLength);
			SAMRecord r = new SAMRecord(getHeader());
			r.setReadName(Integer.toString(i));
			r.setReadBases(B(perfect_read));
			reads.add(r);
		}
		return reads;
	}
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
}
