package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;


public class PathNodeIteratorTest extends TestHelper {
	@Test
	public void should_merge_successive_nodes() {
		List<KmerAggregateNode> in = new ArrayList<KmerAggregateNode>();
		in.add(new KmerAggregateNode(K("GAAA"), 1, 1, 2, true));
		in.add(new KmerAggregateNode(K("AAAT"), 1, 2, 3, true));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 10, 10, 4));
		
		assertEquals(1, out.size());
		assertContains(out, 4, "GAAAT", 1, 2, true, 2);
	}
	public void assertCount(int expected, KmerAggregateNode... in) {
		List<KmerAggregateNode> list = new ArrayList<KmerAggregateNode>();
		for (KmerAggregateNode n : in) list.add(n);
		Collections.sort(list, KmerNode.ByStartPosition);
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(list.iterator(), 10, 10, 4));
		assertEquals(expected, out.size());
	}
	@Test
	public void should_only_merge_fully_matching_nodes() {
		assertCount(1,
			new KmerAggregateNode(K("GAAA"), 1, 1, 2, true),
			new KmerAggregateNode(K("AAAT"), 1, 2, 3, true));
		assertCount(2,
				new KmerAggregateNode(K("GAAA"), 1, 1, 2, true),
				new KmerAggregateNode(K("AAAT"), 1, 2, 3, false));
		assertCount(2,
				new KmerAggregateNode(K("GAAA"), 1, 1, 3, true),
				new KmerAggregateNode(K("AAAT"), 1, 2, 3, true));
		assertCount(2,
				new KmerAggregateNode(K("GAAA"), 1, 1, 2, true),
				new KmerAggregateNode(K("AAAT"), 1, 2, 4, true));
		assertCount(1,
				new KmerAggregateNode(K("GAAA"), 1, 1, 2, true),
				new KmerAggregateNode(K("AAAT"), 2, 2, 3, true));
	}
	@Test
	public void should_handle_self_intersection() {
		assertCount(1, new KmerAggregateNode(K("AAAA"), 1, 1, 5, true));
	}
	@Test
	public void should_split_path_upon_fork() {
		List<KmerAggregateNode> in = new ArrayList<KmerAggregateNode>();
		in.add(new KmerAggregateNode(K("TAAA"), 1, 0, 0, true));
		in.add(new KmerAggregateNode(K("AAAA"), 1, 1, 1, true));
		in.add(new KmerAggregateNode(K("AAAA"), 1, 2, 2, true));
		in.add(new KmerAggregateNode(K("AAAT"), 1, 2, 2, true));
		in.add(new KmerAggregateNode(K("AAAG"), 1, 2, 2, true));
		in.add(new KmerAggregateNode(K("AAAC"), 1, 2, 2, true));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 10, 10, 4));
		assertEquals(5, out.size());
		assertContains(out, 4, "TAAAA", 0, 0, true, 2);
		assertContains(out, 4, "AAAA", 2, 2, true, 1);
		assertContains(out, 4, "AAAC", 2, 2, true, 1);
		assertContains(out, 4, "AAAG", 2, 2, true, 1);
		assertContains(out, 4, "AAAT", 2, 2, true, 1);
	}
	@Test
	public void should_concatenate_to_path() {
		List<KmerAggregateNode> in = new ArrayList<KmerAggregateNode>();
		in.add(new KmerAggregateNode(K("TAAA"), 1, 1, 3, true));
		in.add(new KmerAggregateNode(K("AAAC"), 2, 2, 4, true));
		in.add(new KmerAggregateNode(K("AACT"), 3, 3, 5, true));
		in.add(new KmerAggregateNode(K("ACTT"), 4, 4, 6, true));
		in.add(new KmerAggregateNode(K("CTTG"), 5, 5, 7, true));
		in.add(new KmerAggregateNode(K("TTGG"), 6, 6, 8, true));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 10, 10, 4));
		assertEquals(1, out.size());
		assertContains(out, 4, "TAAACTTGG", 1, 3, true, 1+2+3+4+5+6);
	}
	@Test
	public void should_limit_path_length() {
		List<KmerAggregateNode> in = new ArrayList<KmerAggregateNode>();
		in.add(new KmerAggregateNode(K("TAAA"), 1, 1, 3, true));
		in.add(new KmerAggregateNode(K("AAAC"), 2, 2, 4, true));
		in.add(new KmerAggregateNode(K("AACT"), 3, 3, 5, true));
		in.add(new KmerAggregateNode(K("ACTT"), 4, 4, 6, true));
		in.add(new KmerAggregateNode(K("CTTG"), 5, 5, 7, true));
		in.add(new KmerAggregateNode(K("TTGG"), 6, 6, 8, true));
		List<KmerPathNode> out = Lists.newArrayList(new PathNodeIterator(in.iterator(), 10, 3, 4));
		assertEquals(2, out.size());
		assertContains(out, 4, "TAAACT", 1, 3, true, 1+2+3);
		assertContains(out, 4, "ACTTGG", 4, 6, true, 4+5+6);
	}
	private void assertContains(List<KmerPathNode> list, int k, String baseSequence, int firstStart, int firstEnd, boolean reference, int weight) {
		for (KmerPathNode node : list) {
			if (matches(node, k, baseSequence, firstStart, firstEnd, reference, weight)) {
				return;
			}
		}
		fail(String.format("Missing (%s, %d, %d, %s, %d)", baseSequence, firstStart, firstEnd, reference, weight));
	}
	private boolean matches(KmerPathNode node, int k, String baseSequence, int firstStart, int firstEnd, boolean reference, int weight) {
		if (firstStart != node.startPosition(0)) return false;
		if (firstEnd != node.endPosition(0)) return false;
		if (weight != node.weight()) return false;
		if (node.isReference() != reference) return false;
		List<Long> list = toKmers(baseSequence, k);
		if (list.size() != node.length()) return false;
		for (int i = 0; i < list.size(); i++) {
			if ((long)list.get(i) != node.kmer(i)) return false;
		}
		return true;
	}
	private List<Long> toKmers(String baseSequence, int k) {
		List<Long> list = new ArrayList<Long>();
		for (int i = 0; i + k <= baseSequence.length(); i++) {
			list.add(K(baseSequence.substring(i, i + k)));
		}
		return list;
	}
}
