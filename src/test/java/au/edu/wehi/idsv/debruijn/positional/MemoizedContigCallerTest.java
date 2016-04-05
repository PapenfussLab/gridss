package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.api.client.util.Lists;


public class MemoizedContigCallerTest extends ContigCallerTest {
	@Override
	public ContigCaller getCaller(Iterable<KmerPathNode> input, int maxEvidenceWidth) {
		ArrayList<KmerPathNode> x = Lists.newArrayList(input);
		x.sort(KmerNodeUtil.ByFirstStart);
		return new MemoizedContigCaller(x.iterator(), maxEvidenceWidth);
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
		MemoizedContigCaller caller = new MemoizedContigCaller(input.iterator(), 10);
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
	public void add_should_recalculate_self() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 1, false));
		input.add(KPN(k, "AAAA", 2, 2, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		MemoizedContigCaller caller = new MemoizedContigCaller(input.iterator(), 10);
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
		MemoizedContigCaller caller = new MemoizedContigCaller(input.iterator(), 10);
		assertEquals(2, caller.bestContig().size());
		KmerPathNode added = KPN(k, "AAAA", 3, 3, false);
		KmerPathNode.addEdge(input.get(1), added);
		caller.add(added);
		assertEquals(3, caller.bestContig().size());
	}
}
