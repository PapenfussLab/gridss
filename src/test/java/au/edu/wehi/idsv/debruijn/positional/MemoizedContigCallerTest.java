package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import com.google.api.client.util.Lists;
import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;


public class MemoizedContigCallerTest extends ContigCallerTest {
	@Override
	public ContigCaller getCaller(Iterable<KmerPathNode> input, int maxEvidenceWidth) {
		ArrayList<KmerPathNode> x = Lists.newArrayList(input);
		x.sort(KmerNodeUtil.ByFirstStart);
		return new MemoizedContigCaller(Iterators.peekingIterator(x.iterator()), maxEvidenceWidth);
	}
	@Test
	public void remove_should_recalculate_descendants() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 1, false));
		input.add(KPN(k, "AAAA", 2, 2, false));
		input.add(KPN(k, "AAAA", 3, 3, false));
		input.add(KPN(k, "TTTT", 3, 3, true));
		input.add(KPN(k, "AAAA", 4, 4, false));
		input.add(KPN(k, "AAAA", 5, 5, false));
		input.add(KPN(k, "AAAA", 6, 6, false));
		input.add(KPN(k, "AAAA", 7, 7, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(1), input.get(3));
		KmerPathNode.addEdge(input.get(2), input.get(4));
		KmerPathNode.addEdge(input.get(3), input.get(4));
		KmerPathNode.addEdge(input.get(4), input.get(5));
		KmerPathNode.addEdge(input.get(5), input.get(6));
		KmerPathNode.addEdge(input.get(6), input.get(7));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		ArrayDeque<KmerPathSubnode> pre = caller.bestContig();
		// length 3 with ref anchor
		assertEquals(4, pre.size());
		KmerPathNode ref = input.get(3);
		caller.remove(ref);
		ref.remove();
		ArrayDeque<KmerPathSubnode> post = caller.bestContig();
		assertEquals(7, post.size());
	}
	@Test
	public void remove_should_not_recalculate_self_descendant() {
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(4, "AAAA", 1, 10, false));
		KmerPathNode.addEdge(input.get(0), input.get(0));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		caller.bestContig();
		caller.sanityCheck();
		caller.remove(input.get(0));
		caller.sanityCheck();
		assertNull(caller.bestContig());
	}
	@Test
	public void remove_should_not_add_removed_descendants() {
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(4, "AAAA", 1, 1, false));
		input.add(KPN(4, "AAAA", 2, 2, false));
		input.add(KPN(4, "AAAA", 3, 3, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		caller.bestContig();
		caller.remove(input.get(1));
		caller.remove(input.get(0));
		caller.remove(input.get(2));
		caller.sanityCheck();
		assertNull(caller.bestContig());
	}
	@Test
	public void setremove_should_recalculate_descendants() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 1, false));
		input.add(KPN(k, "AAAA", 2, 2, false));
		input.add(KPN(k, "AAAA", 3, 3, false));
		input.add(KPN(k, "TTTT", 3, 3, true));
		input.add(KPN(k, "AAAA", 4, 4, false));
		input.add(KPN(k, "AAAA", 5, 5, false));
		input.add(KPN(k, "AAAA", 6, 6, false));
		input.add(KPN(k, "AAAA", 7, 7, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(1), input.get(3));
		KmerPathNode.addEdge(input.get(2), input.get(4));
		KmerPathNode.addEdge(input.get(3), input.get(4));
		KmerPathNode.addEdge(input.get(4), input.get(5));
		KmerPathNode.addEdge(input.get(5), input.get(6));
		KmerPathNode.addEdge(input.get(6), input.get(7));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		ArrayDeque<KmerPathSubnode> pre = caller.bestContig();
		// length 3 with ref anchor
		assertEquals(4, pre.size());
		KmerPathNode ref = input.get(3);
		caller.remove(Sets.newHashSet(ref));
		ref.remove();
		ArrayDeque<KmerPathSubnode> post = caller.bestContig();
		assertEquals(7, post.size());
	}
	@Test
	public void setremove_should_not_recalculate_self_descendant() {
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(4, "AAAA", 1, 10, false));
		KmerPathNode.addEdge(input.get(0), input.get(0));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		caller.bestContig();
		caller.sanityCheck();
		caller.remove(Sets.newHashSet(input.get(0)));
		caller.sanityCheck();
		assertNull(caller.bestContig());
	}
	@Test
	public void setremove_should_not_add_removed_descendants() {
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(4, "AAAA", 1, 1, false));
		input.add(KPN(4, "AAAA", 2, 2, false));
		input.add(KPN(4, "AAAA", 3, 3, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		caller.bestContig();
		caller.remove(Sets.newHashSet(input));
		caller.sanityCheck();
		assertNull(caller.bestContig());
	}
	@Test
	public void setremove_should_restart_alternate_paths() {
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(4, "AAAA", 1, 1, true, 2));
		input.add(KPN(4, "AAAA", 2, 2, false));
		input.add(KPN(4, "TTTT", 2, 2, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		caller.bestContig();
		assertEquals(0, caller.tracking_frontierSize());
		caller.remove(Sets.newHashSet(input.get(0), input.get(1)));
		Assert.assertNotNull(caller.bestContig());
	}
	@Test
	public void add_should_recalculate_self() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 1, false));
		input.add(KPN(k, "AAAA", 2, 2, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		assertEquals(2, caller.bestContig().size());
		KmerPathNode added = KPN(k, "AAAA", 0, 0, false);
		KmerPathNode.addEdge(added, input.get(0));
		caller.add(added);
		assertEquals(3, caller.bestContig().size());
	}
	@Test
	public void add_should_recalculate_pre() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 1, false));
		input.add(KPN(k, "AAAA", 2, 2, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(input.iterator()), 10);
		assertEquals(2, caller.bestContig().size());
		KmerPathNode added = KPN(k, "AAAA", 3, 3, false);
		KmerPathNode.addEdge(input.get(1), added);
		caller.add(added);
		assertEquals(3, caller.bestContig().size());
	}
	/**
	 * Turns out the frontier parent is not necessarily memoized.
	 */
	@Test
	public void not_invariant_frontier_parent_object_is_memoized() {
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(new ArrayList<KmerPathNode>().iterator()), 10);
		KmerPathNode n1 = KPN(1, "A", 1, 2, false);
		KmerPathNode n2 = KPN(1, "C", 2, 3, false);
		KmerPathNode n3 = KPN(1, "G", 3, 3, false);
		KmerPathNode.addEdge(n1, n2);
		KmerPathNode.addEdge(n2, n3);
		caller.add(n1);
		caller.add(n2);
		caller.add(n3);
		assertEquals(3, caller.bestContig().size());
		// now we break the n2 best path into new best paths
		KmerPathNode n4 = KPN(1, "T", 2, 2, false, 10);
		KmerPathNode.addEdge(n4, n2);
		caller.add(n4);
		assertEquals(n4, caller.bestContig().getFirst().node());
		// check that set removed the memoized path at n3
		caller.remove(Sets.newHashSet(n1, n4));
		assertEquals(2, caller.bestContig().size());
	}
	@Test
	public void remove_should_remove_all_descendents() {
		MemoizedContigCaller caller = new MemoizedContigCaller(Iterators.peekingIterator(new ArrayList<KmerPathNode>().iterator()), 10);
		KmerPathNode n1 = KPN(1, "A", 1, 2, false);
		KmerPathNode n2 = KPN(1, "C", 2, 3, false);
		KmerPathNode n3 = KPN(1, "G", 3, 3, false);
		KmerPathNode.addEdge(n1, n2);
		KmerPathNode.addEdge(n2, n3);
		caller.add(n1);
		caller.add(n2);
		caller.add(n3);
		assertEquals(3, caller.bestContig().size());
		// now we break the n2 best path into new best paths
		KmerPathNode n4 = KPN(1, "T", 2, 2, false, 10);
		KmerPathNode.addEdge(n4, n2);
		caller.add(n4);
		assertEquals(n4, caller.bestContig().getFirst().node());
		// testing 
		caller.remove(n1);
		caller.remove(n4);
		// check that set removed the memoized path at n3
		assertEquals(2, caller.bestContig().size());
	}
}
