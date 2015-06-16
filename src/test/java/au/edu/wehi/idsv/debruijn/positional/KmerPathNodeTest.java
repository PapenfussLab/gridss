package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import org.junit.Assert;
import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class KmerPathNodeTest extends TestHelper {
	@Test
	public void Constructor_should_copy_KmerNode() {
		KmerPathNode pn = new KmerPathNode(new ImmutableKmerNode(0, 1, 2, 3, true));
		assertEquals(1, pn.startPosition());
		assertEquals(2, pn.endPosition());
		assertEquals(0, pn.kmer());
		assertEquals(3, pn.weight());
		assertTrue(pn.isReference());
		assertTrue(pn.isValid());
	}
	@Test
	public void append_should_add_to_end() {
		ImmutableKmerNode n1 = new ImmutableKmerNode(0, 1, 2, 2, false);
		ImmutableKmerNode n2 = new ImmutableKmerNode(1, 2, 3, 3, false);
		ImmutableKmerNode n3 = new ImmutableKmerNode(2, 3, 4, 4, false);
		
		KmerPathNode pn = new KmerPathNode(n1);
		pn.append(n2);
		pn.append(n3);
		assertEquals(2+3+4, pn.weight());
		assertFalse(pn.isReference());
		assertTrue(pn.isValid());
		
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
		ImmutableKmerNode n0 = new ImmutableKmerNode(0, 1, 2, 2, false);
		ImmutableKmerNode n1 = new ImmutableKmerNode(0, 1, 2, 2, false);
		ImmutableKmerNode n2 = new ImmutableKmerNode(1, 2, 3, 3, false);
		ImmutableKmerNode n3 = new ImmutableKmerNode(2, 3, 4, 4, false);
		ImmutableKmerNode n4 = new ImmutableKmerNode(3, 4, 5, 6, false);
		ImmutableKmerNode n5 = new ImmutableKmerNode(0, 1, 2, 2, false);
		
		KmerPathNode pn0 = new KmerPathNode(n0);
		KmerPathNode pn0a = new KmerPathNode(new ImmutableKmerNode(-1, 1, 2, 2, false));
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
	public KmerPathNode kpn(long[] kmers, int[] weights, int start, int end, boolean reference) {
		KmerPathNode pn = new KmerPathNode(kmers[0], weights[0], start, end, reference);
		for (int i = 1; i < kmers.length; i++) {
			pn.append(new ImmutableKmerNode(kmers[i], start + i, end + i, weights[i], reference));
		}
		return pn;
	}
	@Test
	public void canCoalese_should_require_adjacent_before_position_and_everything_else_matching() {
		assertTrue(kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 5, 10, true).canCoalese(
				   kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 4, true)));
		
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 5, 10, true).canCoalese(
					kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 3, true)));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 3, 3, 4 }, 5, 10, true).canCoalese(
				    kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 4, true)));
		assertFalse(kpn(new long[] { 0, 3, 2, 3 }, new int[] { 1, 2, 3, 4 }, 5, 10, true).canCoalese(
				    kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 4, true)));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 5, 10, true).canCoalese(
				    kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 4, false)));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 4, true).canCoalese(
				    kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 5, 10, true)));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 5, 10, true).canCoalese(
				    kpn(new long[] { 0, 1, 2, 3, 4 }, new int[] { 1, 2, 3, 4, 5 }, 3, 4, true)));
		assertFalse(kpn(new long[] { 0, 1, 2, 3 }, new int[] { 4, 3, 2, 1 }, 5, 10, true).canCoalese(
				    kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 4, true)));
	}
	@Test
	public void coaleseAdjacent_should_merge_interval() {
		KmerPathNode pn = kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 5, 10, true);
		KmerPathNode pn2 = kpn(new long[] { 0, 1, 2, 3 }, new int[] { 1, 2, 3, 4 }, 3, 4, true);
		pn.coaleseAdjacent(pn2);
		assertEquals(3, pn.startPosition(0));
		assertEquals(10, pn.endPosition(0));
		assertEquals(0, pn.kmer(0));
		assertEquals(1+2+3+4, pn.weight());
		assertTrue(pn.isReference());
		assertTrue(pn.isValid());
		assertFalse(pn2.isValid());
	}
	@Test
	public void KmerNode_kmer_position_should_be_of_lastKmer() {
		ImmutableKmerNode n1 = new ImmutableKmerNode(0, 1, 2, 2, false);
		ImmutableKmerNode n2 = new ImmutableKmerNode(1, 2, 3, 3, false);
		ImmutableKmerNode n3 = new ImmutableKmerNode(2, 3, 4, 4, false);
		
		KmerPathNode pn = new KmerPathNode(n1);
		pn.append(n2);
		pn.append(n3);
		assertEquals(pn.startPosition(2), pn.startPosition());
		assertEquals(pn.endPosition(2), pn.endPosition());
		assertEquals(pn.kmer(2), pn.kmer());
	}
	@Test
	public void coaleseAdjacent_should_merge_edges() {
		// TODO: Assert.fail();
	}	
}
