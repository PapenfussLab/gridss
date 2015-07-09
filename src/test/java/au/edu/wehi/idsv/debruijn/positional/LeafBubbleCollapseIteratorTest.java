package au.edu.wehi.idsv.debruijn.positional;

import java.util.Iterator;


public class LeafBubbleCollapseIteratorTest extends CollapseIteratorTest {
	@Override
	protected CollapseIterator create(Iterator<KmerPathNode> it, int k,
			int maxPathCollapseLength, int maxBasesMismatch,
			boolean bubblesAndLeavesOnly, double minimumPathNodeEntropy) {
		return new LeafBubbleCollapseIterator(it, k, maxPathCollapseLength, maxBasesMismatch);
	}
}
