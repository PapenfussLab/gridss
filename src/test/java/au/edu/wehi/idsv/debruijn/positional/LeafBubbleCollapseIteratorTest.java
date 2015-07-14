package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;


public class LeafBubbleCollapseIteratorTest extends CollapseIteratorTest {
	@Override
	protected CollapseIterator create(Iterator<KmerPathNode> it, int k,
			int maxPathCollapseLength, int maxBasesMismatch,
			boolean bubblesAndLeavesOnly, double minimumPathNodeEntropy) {
		return new LeafBubbleCollapseIterator(it, k, maxPathCollapseLength, maxBasesMismatch);
	}
	@Test
	public void should_not_merge_non_bubble() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAT", 1, 1, false)); 
		input.add(KPN(k, "AATC", 2, 2, false));
		input.add(KPN(k,  "ATCCAAA", 3, 3, false));
		input.add(KPN(k, "AATGGAAA", 2, 2, false));
		input.add(KPN(k, "AAACTTTTTTTT", 7, 7, false));
		input.add(KPN(k, "AAAGGGGGGGGG", 7, 7, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(0), input.get(3));
		KmerPathNode.addEdge(input.get(2), input.get(4));
		KmerPathNode.addEdge(input.get(2), input.get(5));
		KmerPathNode.addEdge(input.get(3), input.get(4));
		KmerPathNode.addEdge(input.get(3), input.get(5));
		List<KmerPathNode> result = go(k, 100, 2, input);
		assertEquals(6, result.size());
	}
}
