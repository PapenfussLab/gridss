package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;

import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;


public class PathCollapseIteratorTest extends CollapseIteratorTest {
	@Override
	protected CollapseIterator create(Iterator<KmerPathNode> it, int k,
			int maxPathCollapseLength, int maxBasesMismatch,
			boolean bubblesAndLeavesOnly, double minimumPathNodeEntropy) {
		// TODO Auto-generated method stub
		return new PathCollapseIterator(it, k, maxPathCollapseLength, maxBasesMismatch, bubblesAndLeavesOnly, minimumPathNodeEntropy);
	}
	@Test
	public void ISSUE_merge_collapse_non_bubble_results_in_invalid_kmer_path() {
		KmerPathNode start = KPN(4, "AAAA", 1, 1, true);
		KmerPathNode toMerge1 = KPN(4, "AAAT", 2, 2, true);
		KmerPathNode toMerge2 = KPN(4, "AATGGG", 3, 3, true);
		KmerPathNode into = KPN(4, "AAACGGG", 2, 2, true);
		KmerPathNode end = KPN(4, "GGGG", 6, 6, true);
		KmerPathNode toMergeAlt = KPN(4, "AATCCCGTGTGTTGTG", 3, 3, true);
		KmerPathNode.addEdge(start, toMerge1);
		KmerPathNode.addEdge(toMerge1, toMerge2);
		KmerPathNode.addEdge(start, into);
		KmerPathNode.addEdge(toMerge2, end);
		KmerPathNode.addEdge(into, end);
		KmerPathNode.addEdge(toMerge1, toMergeAlt);
		List<KmerPathNode> result = go(4, 100, 2, ImmutableList.of(start, toMerge1, toMerge2, into, end, toMergeAlt));
		assertEquals(5, result.size());
		// problem:
		// start > into(toMerge1) > toMergeAlt path kmers are inconsistent!
		// AAAATCCCGTGTGTTGTG
		//     |
		// AAAACCCCGTGTGTTGTG
	}
}
