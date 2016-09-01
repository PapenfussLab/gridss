package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class MemoizedTraverseTest extends TestHelper {
	public KmerPathNode N(int start, int stop, int weight) {
		return KPN(1, "A", start, stop, false, weight);
	}
	@Test
	public void should_memoize_best_path() {
		KmerPathNode n = N(1, 1, 1);
		TraversalNode tn = new TraversalNode(new KmerPathSubnode(n), 0);
		MemoizedTraverse mt = new MemoizedTraverse();
		mt.memoize(tn);
		assertTrue(mt.memoized(n).contains(tn));
	}
	@Test
	public void should_not_memoize_suboptimal_path() {
		KmerPathNode n = N(1, 1, 1);
		TraversalNode tn = new TraversalNode(new KmerPathSubnode(n), 1);
		TraversalNode tn2 = new TraversalNode(new KmerPathSubnode(n), 0);
		MemoizedTraverse mt = new MemoizedTraverse();
		mt.memoize(tn);
		mt.memoize(tn2);
		assertTrue(mt.memoized(n).contains(tn));
		assertFalse(mt.memoized(n).contains(tn2));
	}
	@Test
	public void should_add_memoized_path_to_frontier() {
		KmerPathNode n = N(1, 1, 1);
		TraversalNode tn = new TraversalNode(new KmerPathSubnode(n), 0);
		MemoizedTraverse mt = new MemoizedTraverse();
		mt.memoize(tn);
		assertEquals(tn, mt.peekFrontier());
	}
	@Test
	public void new_best_path_should_split_old() {
		KmerPathNode n = N(1, 20, 1);
		TraversalNode tn = new TraversalNode(new KmerPathSubnode(n, 1, 8), 0);
		TraversalNode tn2 = new TraversalNode(new KmerPathSubnode(n, 3, 4), 1);
		MemoizedTraverse mt = new MemoizedTraverse();
		mt.memoize(tn);
		mt.memoize(tn2);
		ArrayList<TraversalNode> result = Lists.newArrayList(mt.memoized(n));
		assertEquals(3, result.size());
		assertNodeIs(1, 2, 1, result.get(0));
		assertNodeIs(3, 4, 2, result.get(1));
		assertNodeIs(5, 8, 1, result.get(2));
		assertEquals(1, mt.pollFrontier().node.firstStart());
		assertEquals(3, mt.pollFrontier().node.firstStart());
		assertEquals(5, mt.pollFrontier().node.firstStart());
		assertTrue(mt.isEmptyFrontier());
	}
	private void assertNodeIs(int start, int end, int score, TraversalNode n) {
		assertEquals(start, n.node.firstStart());
		assertEquals(end, n.node.firstEnd());
		assertEquals(score, n.score);
	}
	@Test
	public void new_path_should_be_split_if_not_best() {
		KmerPathNode n = N(1, 20, 1);
		MemoizedTraverse mt = new MemoizedTraverse();
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 3, 3), 1));
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 4, 4), 2));
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 7, 7), 3));
		mt.pollFrontier();
		mt.pollFrontier();
		mt.pollFrontier();
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 1, 8), 0));
		ArrayList<TraversalNode> result = Lists.newArrayList(mt.memoized(n));
		assertNodeIs(1, 2, 1, result.get(0));
		assertNodeIs(3, 3, 2, result.get(1));
		assertNodeIs(4, 4, 3, result.get(2));
		assertNodeIs(5, 6, 1, result.get(3));
		assertNodeIs(7, 7, 4, result.get(4));
		assertNodeIs(8, 8, 1, result.get(5));
		assertEquals(1, mt.pollFrontier().node.firstStart());
		assertEquals(5, mt.pollFrontier().node.firstStart());
		assertEquals(8, mt.pollFrontier().node.firstStart());
		assertTrue(mt.isEmptyFrontier());
	}
	@Test
	public void new_path_should_be_split_if_not_best_at_edge() {
		KmerPathNode n = N(1, 20, 1);
		MemoizedTraverse mt = new MemoizedTraverse();
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 3, 3), 1));
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 4, 4), 2));
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 7, 7), 3));
		mt.pollFrontier();
		mt.pollFrontier();
		mt.pollFrontier();
		mt.memoize(new TraversalNode(new KmerPathSubnode(n, 3, 7), 0));
		ArrayList<TraversalNode> result = Lists.newArrayList(mt.memoized(n));
		assertNodeIs(3, 3, 2, result.get(0));
		assertNodeIs(4, 4, 3, result.get(1));
		assertNodeIs(5, 6, 1, result.get(2));
		assertNodeIs(7, 7, 4, result.get(3));
		assertEquals(5, mt.pollFrontier().node.firstStart());
		assertTrue(mt.isEmptyFrontier());
	}
	@Test
	public void frontier_should_be_sorted_by_end_position_of_first_kmer() {
		KmerPathNode n1 = N(1, 20, 1);
		KmerPathNode n2 = N(1, 20, 1);
		MemoizedTraverse mt = new MemoizedTraverse();
		mt.memoize(new TraversalNode(new KmerPathSubnode(n1, 1, 10), 0));
		mt.memoize(new TraversalNode(new KmerPathSubnode(n2, 2, 5), 0));
		assertEquals(2, mt.pollFrontier().node.firstStart());
		assertEquals(1, mt.pollFrontier().node.firstStart());
		assertTrue(mt.isEmptyFrontier());
	}
	@Test
	public void unmemoize_should_not_remove_path_if_outdated_parent_overlaps_removal_interval() {
		// 1234
		// a
		// \
		//  bb
		//   \ 
		//    c
		// 1234
		KmerPathNode a = N(1, 1, 2);
		KmerPathNode b = N(2, 3, 1);
		KmerPathNode c = N(4, 4, 1);
		KmerPathNode.addEdge(a, b);
		KmerPathNode.addEdge(b, c);
		MemoizedTraverse mt = new MemoizedTraverse();
		TraversalNode b1 = new TraversalNode(new KmerPathSubnode(b, 2, 3), 0);
		TraversalNode bc1 = new TraversalNode(b1, new KmerPathSubnode(c, 4, 4));
		mt.memoize(b1);
		mt.memoize(bc1);
		TraversalNode a2 = new TraversalNode(new KmerPathSubnode(a, 1, 1), 0);
		TraversalNode ab2 = new TraversalNode(a2, new KmerPathSubnode(b, 2, 2));
		mt.memoize(a2);
		mt.memoize(ab2);
		while (!mt.isEmptyFrontier()) mt.pollFrontier();
		
		// ok, now we're set up:
		// the memoized value at c points to an old b frontier that
		// overlaps the a path
		mt.remove(a);
		
		// we still shouldn't remove c because the actual path doesn't overlap
		assertEquals(1, mt.memoized(c).size());
	}
	@Test
	public void unmemoize_should_not_remove_path_if_outdated_descendant_overlaps_removal_interval() {
		// 123456
		// a
		// \
		//  bb
		//  \\
		//   cc
		//    \ 
		//     d
		// 12345
		KmerPathNode a = N(1, 2, 2);
		KmerPathNode b = N(2, 3, 1);
		KmerPathNode c = N(3, 4, 1);
		KmerPathNode d = N(5, 5, 1);
		KmerPathNode.addEdge(a, b);
		KmerPathNode.addEdge(b, c);
		KmerPathNode.addEdge(c, d);
		MemoizedTraverse mt = new MemoizedTraverse();
		TraversalNode b1 = new TraversalNode(new KmerPathSubnode(b, 2, 3), 0);
		TraversalNode bc1 = new TraversalNode(b1, new KmerPathSubnode(c, 3, 4));
		TraversalNode bcd1 = new TraversalNode(bc1, new KmerPathSubnode(d, 5, 5));
		mt.memoize(b1);
		mt.memoize(bc1);
		mt.memoize(bcd1);
		TraversalNode a2 = new TraversalNode(new KmerPathSubnode(a, 1, 1), 0);
		TraversalNode ab2 = new TraversalNode(a2, new KmerPathSubnode(b, 2, 2));
		TraversalNode abc2 = new TraversalNode(ab2, new KmerPathSubnode(c, 3, 3));
		mt.memoize(a2);
		mt.memoize(ab2);
		mt.memoize(abc2);
		while (!mt.isEmptyFrontier()) mt.pollFrontier();
		
		// ok, now we're set up:
		// the memoized value at d points to an old c frontier that
		// overlaps the a path
		mt.remove(a);
		
		// we still shouldn't remove c because the actual path doesn't overlap
		assertEquals(1, mt.memoized(d).size());
	}
}
