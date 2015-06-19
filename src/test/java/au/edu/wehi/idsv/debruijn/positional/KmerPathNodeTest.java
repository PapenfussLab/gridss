package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.stream.IntStream;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class KmerPathNodeTest extends TestHelper {
	@Test
	public void Constructor_should_copy_KmerNode() {
		KmerPathNode pn = new KmerPathNode(new ImmutableKmerNode(0, 2, 3, true, 1));
		assertEquals(2, pn.startPosition());
		assertEquals(3, pn.endPosition());
		assertEquals(0, pn.kmer());
		assertEquals(1, pn.weight());
		assertTrue(pn.isReference());
		assertTrue(pn.isValid());
	}
	@Test
	public void append_should_add_to_end() {
		ImmutableKmerNode n1 = new ImmutableKmerNode(0, 1, 2, false, 2);
		ImmutableKmerNode n2 = new ImmutableKmerNode(1, 2, 3, false, 3);
		ImmutableKmerNode n3 = new ImmutableKmerNode(2, 3, 4, false, 4);
		
		KmerPathNode pn = new KmerPathNode(n1);
		pn.append(n2);
		pn.append(n3);
		assertEquals(2+3+4, pn.weight());
		assertFalse(pn.isReference());
		assertTrue(pn.isValid());
		assertEquals(3, pn.length());
		assertEquals(2, pn.width());
		
		assertEquals(1, pn.startPosition(0));
		assertEquals(2, pn.endPosition(0));
		assertEquals(0, pn.kmer(0));
		assertEquals(2, pn.weight(0));
		assertEquals(2, pn.startPosition(1));
		assertEquals(3, pn.endPosition(1));
		assertEquals(1, pn.kmer(1));
		assertEquals(3, pn.weight(1));
		assertEquals(3, pn.startPosition(2));
		assertEquals(4, pn.endPosition(2));
		assertEquals(2, pn.kmer(2));
		assertEquals(4, pn.weight(2));
	}
	@Test
	public void prepend_should_add_to_start() {
		ImmutableKmerNode n0 = new ImmutableKmerNode(0, 1, 2, false, 2);
		ImmutableKmerNode n1 = new ImmutableKmerNode(0, 1, 2, false, 2);
		ImmutableKmerNode n2 = new ImmutableKmerNode(1, 2, 3, false, 3);
		ImmutableKmerNode n3 = new ImmutableKmerNode(2, 3, 4, false, 4);
		ImmutableKmerNode n4 = new ImmutableKmerNode(3, 4, 5, false, 6);
		ImmutableKmerNode n5 = new ImmutableKmerNode(0, 1, 2, false, 2);
		
		KmerPathNode pn0 = new KmerPathNode(n0);
		KmerPathNode pn0a = new KmerPathNode(new ImmutableKmerNode(-1, 2, 2, false, 1));
		KmerPathNode pn5 = new KmerPathNode(n5);
		
		KmerPathNode pn1 = new KmerPathNode(n1);
		pn1.append(n2);
		KmerPathNode pn2 = new KmerPathNode(n3);
		pn2.append(n4);
		KmerPathNode.addEdge(pn0a, pn1);
		KmerPathNode.addEdge(pn0, pn1);
		KmerPathNode.addEdge(pn1, pn2);
		KmerPathNode.addEdge(pn2, pn5);
		
		pn2.prepend(pn1);
		assertFalse(pn1.isValid());
		assertTrue(pn2.next().size() == 1);
		assertTrue(pn2.prev().size() == 2);
		assertEquals(2+3+4+6, pn2.weight());
		assertFalse(pn2.isReference());
		assertTrue(pn2.isValid());
		
		assertEquals(1, pn2.startPosition(0));
		assertEquals(2, pn2.endPosition(0));
		assertEquals(0, pn2.kmer(0));
		assertEquals(2, pn2.weight(0));
		assertEquals(2, pn2.startPosition(1));
		assertEquals(3, pn2.endPosition(1));
		assertEquals(1, pn2.kmer(1));
		assertEquals(3, pn2.weight(1));
		assertEquals(3, pn2.startPosition(2));
		assertEquals(4, pn2.endPosition(2));
		assertEquals(2, pn2.kmer(2));
		assertEquals(4, pn2.weight(2));
		assertEquals(4, pn2.startPosition(3));
		assertEquals(5, pn2.endPosition(3));
		assertEquals(3, pn2.kmer(3));
		assertEquals(6, pn2.weight(3));
	}
	@Test
	public void prepend_should_relink_earlier_nodes() {
		KmerPathNode pn1 = new KmerPathNode(0, 1, 1, true, 1);
		KmerPathNode pn2 = new KmerPathNode(0, 2, 2, true, 1);
		KmerPathNode pn3 = new KmerPathNode(0, 3, 3, true, 1);
		KmerPathNode.addEdge(pn1, pn2);
		KmerPathNode.addEdge(pn2, pn3);
		pn3.prepend(pn2);
		assertEquals(pn3, pn1.next().get(0));
	}
	public KmerPathNode kpn(long[] kmers, int start, int end, boolean reference, int[] weights) {
		KmerPathNode pn = new KmerPathNode(kmers[0], start, end, reference, weights[0]);
		for (int i = 1; i < kmers.length; i++) {
			pn.append(new ImmutableKmerNode(kmers[i], start + i, end + i, reference, weights[i]));
		}
		return pn;
	}
	public KmerPathNode kpn(long[] kmers, int start, int end, boolean reference) {
		int[] weights = new int[kmers.length];
		Arrays.fill(weights, 1);
		return kpn(kmers, start, end, reference, weights);
	}
	public void assertIs(KmerPathNode pn, long[] kmers, int start, int end, boolean reference, int[] weights) {
		assertEquals(start, pn.startPosition(0));
		assertEquals(end, pn.endPosition(0));
		assertEquals(kmers.length, pn.length());
		for (int i = 0; i < pn.length(); i++) {
			assertEquals(kmers[i], pn.kmer(i));
		}
		for (int i = 0; i < pn.length(); i++) {
			assertEquals(weights[i], pn.weight(i));
		}
		assertEquals(IntStream.of(weights).sum(), pn.weight());
		assertEquals(reference, pn.isReference());
		assertTrue(pn.isValid());
		for (int i = 1; i < pn.next().size(); i++) {
			assertTrue(pn.next().get(i - 1).startPosition(0) <= pn.next().get(i).startPosition(0));
		}
		for (int i = 1; i < pn.prev().size(); i++) {
			assertTrue(pn.prev().get(i - 1).startPosition() <= pn.prev().get(i).startPosition());
		}
	}
	@Test
	public void canCoalese_should_require_adjacent_before_position_and_everything_else_matching() {
		assertTrue(kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 }).canCoaleseBeforeAdjacent(
				   kpn(new long[] { 0, 1, 2, 3 }, 3, 4, true, new int[] { 1, 2, 3, 4 })));
		
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 }).canCoaleseBeforeAdjacent(
					kpn(new long[] { 0, 1, 2, 3 }, 3, 3, true, new int[] { 1, 2, 3, 4 })));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 3, 3, 4 }).canCoaleseBeforeAdjacent(
				    kpn(new long[] { 0, 1, 2, 3 }, 3, 4, true, new int[] { 1, 2, 3, 4 })));
		assertFalse(kpn(new long[] { 0, 3, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 }).canCoaleseBeforeAdjacent(
				    kpn(new long[] { 0, 1, 2, 3 }, 3, 4, true, new int[] { 1, 2, 3, 4 })));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 }).canCoaleseBeforeAdjacent(
				    kpn(new long[] { 0, 1, 2, 3 }, 3, 4, false, new int[] { 1, 2, 3, 4 })));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, 3, 4, true, new int[] { 1, 2, 3, 4 }).canCoaleseBeforeAdjacent(
				    kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 })));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 }).canCoaleseBeforeAdjacent(
				    kpn(new long[] { 0, 1, 2, 3, 4 }, 3, 4, true, new int[] { 1, 2, 3, 4, 5 })));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 4, 3, 2, 1 }).canCoaleseBeforeAdjacent(
				    kpn(new long[] { 0, 1, 2, 3 }, 3, 4, true, new int[] { 1, 2, 3, 4 })));
	}
	@Test
	public void coaleseAdjacent_should_merge_interval() {
		KmerPathNode pn1 = kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 });
		KmerPathNode pn2 = kpn(new long[] { 0, 1, 2, 3 }, 3, 4, true, new int[] { 1, 2, 3, 4 });
		
		KmerPathNode pre1 = kpn(new long[] { 0 }, 4, 5, true, new int[] { 1});
		KmerPathNode.addEdge(pre1, pn1);
		KmerPathNode.addEdge(kpn(new long[] { 0 }, 6, 7, true, new int[] { 1}), pn1);
		KmerPathNode sharedPrev = kpn(new long[] { 1 }, 0, 6, true, new int[] { 1});
		KmerPathNode.addEdge(sharedPrev, pn1);
		KmerPathNode.addEdge(sharedPrev, pn2);
		KmerPathNode pre2 = kpn(new long[] { 2 }, 2, 2, true, new int[] { 1 });
		KmerPathNode.addEdge(pre2, pn2);
		KmerPathNode sharedNext = kpn(new long[] { 4 }, 0, 20, true, new int[] { 1 });
		KmerPathNode.addEdge(pn1, sharedNext);
		KmerPathNode.addEdge(pn2, sharedNext);
		
		pn1.coaleseBeforeAdjacent(pn2);
		assertEquals(4, pn1.prev().size());
		assertEquals(1, pn1.next().size());
		assertEquals(1, sharedPrev.next().size());
	}
	@Test
	public void invalidate_should_remove_node() {
		KmerPathNode pn1 = kpn(new long[] { 0, 1, 2, 3 }, 5, 10, true, new int[] { 1, 2, 3, 4 });
		pn1.invalidate();
		assertFalse(pn1.isValid());
	}
	@Test
	public void next_should_sort_by_first_kmer_start() {
		KmerPathNode pn1 = kpn(new long[] { 0 }, 1, 100, true, new int[] { 1 });
		KmerPathNode.addEdge(pn1, kpn(new long[] { 0, 1, 2, 3, 4, 5, 6, 7 }, 3, 3, true));
		KmerPathNode.addEdge(pn1, kpn(new long[] { 0, 1, 3, 4 }, 2, 2, true));
		assertEquals(2, pn1.next().get(0).startPosition(0));
		assertEquals(3, pn1.next().get(1).startPosition(0));
	}
	@Test
	public void prev_should_sort_by_last_kmer_start() {
		KmerPathNode pn1 = kpn(new long[] { 0 }, 1, 100, true, new int[] { 1 });
		KmerPathNode.addEdge(kpn(new long[] { 0, 1, 2, 3, 4, 5, 6, 7 }, 1, 1, true), pn1);
		KmerPathNode.addEdge(kpn(new long[] { 0, 1, 3, 4 }, 2, 2, true), pn1);
		assertEquals(5, pn1.prev().get(0).startPosition());
		assertEquals(8, pn1.prev().get(1).startPosition());
	}
	@Test
	public void KmerNode_kmer_position_should_be_of_lastKmer() {
		ImmutableKmerNode n1 = new ImmutableKmerNode(0, 2, 2, false, 1);
		ImmutableKmerNode n2 = new ImmutableKmerNode(1, 3, 3, false, 2);
		ImmutableKmerNode n3 = new ImmutableKmerNode(2, 4, 4, false, 3);
		
		KmerPathNode pn = new KmerPathNode(n1);
		pn.append(n2);
		pn.append(n3);
		assertEquals(pn.startPosition(2), pn.startPosition());
		assertEquals(pn.endPosition(2), pn.endPosition());
		assertEquals(pn.kmer(2), pn.kmer());
	}
	@Test
	public void splitAtLength_should_break_after_nth_kmer() {
		KmerPathNode pn = kpn(new long[] { 0, 1, 2, 3, 4 }, 1, 10, true, new int[] { 1, 2, 3, 4, 5 });
		KmerPathNode.addEdge(kpn(new long[] { 0, 1, 2, 3, 4, 5, 6, 7 }, 1, 1, true), pn);
		KmerPathNode.addEdge(kpn(new long[] { 0, 1, 2, 3, 4, 5, 6, 7 }, 2, 2, true), pn);
		KmerPathNode.addEdge(pn, kpn(new long[] { 0, 1, 2, 3, 4, 5, 6, 7 }, 1, 1, false));
		KmerPathNode.addEdge(pn, kpn(new long[] { 0, 1, 2, 3, 4, 5, 6, 7 }, 2, 2, false));
		KmerPathNode split = pn.splitAtLength(3);
		
		assertIs(split, new long[] { 0, 1, 2 }, 1, 10, true, new int[] { 1, 2, 3 });
		assertIs(pn, new long[] { 3, 4 }, 4, 13, true, new int[] { 4, 5 });
		
		assertEquals(1, pn.prev().size());
		assertEquals(1, split.next().size());
		assertEquals(split, pn.prev().get(0));
		assertEquals(pn, split.next().get(0));
		assertEquals(2, split.prev().size());
		assertEquals(2, pn.next().size());
	}
	@Test
	public void splitAtStartPosition_should_split_so_first_kmer_starts_at_given_position() {
		KmerPathNode pn = kpn(new long[] { 0, 1, 2, 3, 4 }, 1, 10, true, new int[] { 1, 2, 3, 4, 5 });
		// prev: 
		KmerPathNode.addEdge(kpn(new long[] { 0, 1 }, -5, 1, true), pn); // split only
		KmerPathNode.addEdge(kpn(new long[] { 1, 1 }, 1, 1, true), pn); // split only
		KmerPathNode.addEdge(kpn(new long[] { 2, 1 }, 2, 2, true), pn); // pn only
		KmerPathNode.addEdge(kpn(new long[] { 3, 1 }, -100, 100, false), pn); // both
		// next
		KmerPathNode.addEdge(pn, kpn(new long[] { 0, 1, 2, 3, 4, 5, 6, 7 }, 6, 6, false)); // split only
		KmerPathNode.addEdge(pn, kpn(new long[] { 1, 1, 2, 3, 4, 5, 6, 7 }, 8, 8, false)); // split only
		KmerPathNode.addEdge(pn, kpn(new long[] { 2, 1, 2, 3, 4, 5, 6, 7 }, 8, 9, false)); // both
		KmerPathNode.addEdge(pn, kpn(new long[] { 3, 1, 2, 3, 4, 5, 6, 7 }, 9, 9, false)); // pn only
		KmerPathNode split = pn.splitAtStartPosition(4);
		
		assertIs(split, new long[] { 0, 1, 2, 3, 4 }, 1, 3, true, new int[] { 1, 2, 3, 4, 5 });
		assertIs(pn, new long[] { 0, 1, 2, 3, 4 }, 4, 10, true, new int[] { 1, 2, 3, 4, 5 });
		
		assertEquals(2, pn.next().size());
		assertEquals(3, split.next().size());
		
		assertEquals(2, pn.prev().size());
		assertEquals(3, split.prev().size());
	}
}
