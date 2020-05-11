package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class AggregateNodeIteratorTest extends TestHelper {
	public void assertIs(KmerNode node, long kmer, int start, int end, int weight, boolean isReference) {
		assertEquals(kmer, node.lastKmer());
		assertEquals(start, node.lastStart());
		assertEquals(end, node.lastEnd());
		assertEquals(weight, node.weight());
		assertEquals(isReference, node.isReference());
	}
	@Test
	public void single_kmer_should_aggregate_nonreference_into_reference() {
		List<KmerNode> input = new ArrayList<KmerNode>();
		input.add(new ImmutableKmerNode(0, 1, 1, true, 3));
		input.add(new ImmutableKmerNode(0, 1, 1, false, 3));
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(1, anList.size());
		assertIs(anList.get(0), 0, 1, 1, 6, true);
	}
	@Test
	public void single_kmer_should_not_merge_nonreference_interval_with_reference() {
		List<KmerNode> input = new ArrayList<KmerNode>();
		input.add(new ImmutableKmerNode(0, 1, 1, true, 3));
		input.add(new ImmutableKmerNode(0, 2, 2, false, 3));
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(2, anList.size());
	}
	/**
	 * Not critical to merge adjacent nodes as we reduce
	 * after simplification
	 */
	@Test
	public void adjacent_same_weight_need_not_merge() {
		List<KmerNode> input = new ArrayList<KmerNode>();
		input.add(new ImmutableKmerNode(0, 1, 1, true, 3));
		input.add(new ImmutableKmerNode(0, 2, 2, true, 3));
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(2, anList.size());
		assertIs(anList.get(0), 0, 1, 1, 3, true);
		assertIs(anList.get(1), 0, 2, 2, 3, true);
	}
	@Test
	public void adjacent_different_weight_should_not_merge() {
		List<KmerNode> input = new ArrayList<KmerNode>();
		input.add(new ImmutableKmerNode(0, 1, 1, true, 1));
		input.add(new ImmutableKmerNode(0, 2, 2, true, 2));
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(2, anList.size());
		assertIs(anList.get(0), 0, 1, 1, 1, true);
		assertIs(anList.get(1), 0, 2, 2, 2, true);
	}
	@Test
	public void should_merge_overlap_contained() {
		//          1         2         3         4
		// 1234567890123456789012345678901234567890
		// -------
		//  -----
		//   ---
		List<KmerNode> input = new ArrayList<KmerNode>();
		input.add(new ImmutableKmerNode(0, 1, 5, true, 2));
		input.add(new ImmutableKmerNode(0, 2, 4, true, 2));
		input.add(new ImmutableKmerNode(0, 3, 3, true, 2));
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(5, anList.size());
		assertIs(anList.get(0), 0, 1, 1, 2, true);
		assertIs(anList.get(1), 0, 2, 2, 4, true);
		assertIs(anList.get(2), 0, 3, 3, 6, true);
		assertIs(anList.get(3), 0, 4, 4, 4, true);
		assertIs(anList.get(4), 0, 5, 5, 2, true);
	}
	@Test
	public void should_merge_overlap() {
		//          1         2         3         4
		// 1234567890123456789012345678901234567890
		// ----
		//  ----
		//   ----
		// ABCCDE
		List<KmerNode> input = new ArrayList<KmerNode>();
		input.add(new ImmutableKmerNode(0, 1, 4, true, 3));
		input.add(new ImmutableKmerNode(0, 2, 5, true, 4));
		input.add(new ImmutableKmerNode(0, 3, 6, true, 5));
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(5, anList.size());
		assertIs(anList.get(0), 0, 1, 1, 3, true);
		assertIs(anList.get(1), 0, 2, 2, 7, true);
		assertIs(anList.get(2), 0, 3, 4, 12, true);
		assertIs(anList.get(3), 0, 5, 5, 9, true);
		assertIs(anList.get(4), 0, 6, 6, 5, true);
	}
	@Test
	public void should_merge_successive_overlap() {
		//          1         2         3         4
		// 1234567890123456789012345678901234567890
		// ----------
		//  ----------
		//          -----------
		List<KmerNode> input = new ArrayList<KmerNode>();
		input.add(new ImmutableKmerNode(0, 1, 10, true, 1));
		input.add(new ImmutableKmerNode(0, 2, 11, true, 2));
		input.add(new ImmutableKmerNode(0, 10, 20, true, 3));
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(5, anList.size());
		assertIs(anList.get(0), 0, 1, 1, 1, true);
		assertIs(anList.get(1), 0, 2, 9, 3, true);
		assertIs(anList.get(2), 0, 10, 10, 6, true);
		assertIs(anList.get(3), 0, 11, 11, 5, true);
		assertIs(anList.get(4), 0, 12, 20, 3, true);
	}
	@Test
	public void should_merge_sequential_overlap() {
		List<KmerNode> input = new ArrayList<KmerNode>();
		for (int i = 0; i < 64; i++) {
			input.add(new ImmutableKmerNode(0, i, i+4, true, 1));
		}
		Collections.sort(input, KmerNodeUtil.ByLastStart);
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(input.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(68, anList.size());
	}
	@Test
	public void should_aggregate_kmer_weights() {
		int k = 4;
		SAMEvidenceSource ses = SES(30, 60);
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(NRRP(ses, withQual(new byte[] { 1,1,1,1,1}, withSequence("ACGTT", DP(0, 1, "5M", false, 1, 1, "5M", false)))));
		input.add(NRRP(ses, withQual(new byte[] { 1,1,1,1,1}, withSequence("ACGTT", DP(0, 2, "5M", false, 1, 1, "5M", false)))));
		Collections.sort(input, DirectedEvidence.ByStartEnd);
		List<KmerSupportNode> snList = Lists.newArrayList(new SupportNodeIterator(k, input.iterator(), ses.getMaxConcordantFragmentSize(), null, false, 0));
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(snList.iterator()));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(anList));
		assertEquals(4, snList.size());
		assertEquals(6, anList.size());
	}
	@Test
	public void total_weight_should_remain_constant() {
		List<KmerSupportNode> snList = Lists.newArrayList(new SupportNodeIterator(4, SupportNodeIteratorTest.scrp(4, "ACGTTATACCG", 30, 60).iterator(), 60, null, false, 0));
		assertTrue(KmerNodeUtil.ByLastStart.isOrdered(snList));
		List<KmerNode> anList = Lists.newArrayList(new AggregateNodeIterator(snList.iterator()));
		assertEquals(totalWeight(snList), totalWeight(anList));
		//		snList.stream().mapToInt(n -> (n.endPosition() - n.startPosition() + 1) * n.weight()).sum(),
		//		anList.stream().mapToInt(n -> (n.endPosition() - n.startPosition() + 1) * n.weight()).sum());
	}
	//@Test // expensive test to run
	public void should_stream_input() {
		AggregateNodeIterator agIt = new AggregateNodeIterator(new SupportNodeIterator(25, new RandomSoftClipIterator(), 100, null, false, 0));
		while (agIt.hasNext()) {
			agIt.next();
		}
	}
}
