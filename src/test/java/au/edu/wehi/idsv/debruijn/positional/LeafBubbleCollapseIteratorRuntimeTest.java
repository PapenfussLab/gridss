package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class LeafBubbleCollapseIteratorRuntimeTest extends TestHelper {
	protected List<KmerPathNode> go(List<DirectedEvidence> input, int k, int maxLength, int maxPathCollapseLength, int maxBasesMismatch) {
		List<KmerPathNode> pnList = CollapseIteratorTest.asCheckedKPN(k, maxLength, input.toArray(new DirectedEvidence[0]));
		return go(k, maxPathCollapseLength, maxBasesMismatch, pnList);
	}
	protected List<KmerPathNode> go(int k, int maxPathCollapseLength, int maxBasesMismatch, List<KmerPathNode> pnList) {
		assertDisjoint(pnList); // precondition: is the test case well formed
		int pnTotalWeight = totalWeight(pnList);
		CollapseIterator pcit = create(pnList.iterator(), k, maxPathCollapseLength, maxBasesMismatch, false, 0);
		ArrayList<KmerPathNode> result = Lists.newArrayList(pcit);
		assertEquals(pnTotalWeight, totalWeight(result));
		return result;
	}
	protected CollapseIterator create(Iterator<KmerPathNode> it, int k,
			int maxPathCollapseLength, int maxBasesMismatch,
			boolean bubblesAndLeavesOnly, double minimumPathNodeEntropy) {
		return new TestLeafBubbleCollapseIterator(it, k, maxPathCollapseLength, maxBasesMismatch);
	}
	private class TestLeafBubbleCollapseIterator extends LeafBubbleCollapseIterator {
		public TestLeafBubbleCollapseIterator(Iterator<KmerPathNode> it, int k, int maxPathCollapseLength, int maxBasesMismatch) {
			super(it, k, maxPathCollapseLength, maxBasesMismatch);
		}
		@Override
		public boolean collapse(KmerPathNode node, int maxCollapseLength) {
			if (node.firstKmer() != K("AAAT")) return false;
			return super.collapse(node, maxCollapseLength);
		}
	}
	@Test(timeout=2000)
	public void should_not_take_exponential_time() {
		int k = 4;
		int loops = 128;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 1, false)); // 0
		input.add(KPN(k, "AAATTTT", 2, 2, false)); // 1
		input.add(KPN(k, "AAAGAAA", 2, 2, false)); // 2
		input.add(KPN(k, "AAACAAA", 2, 2, false)); // 3
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(0), input.get(3));
		for (int i = 0; i < loops; i++) {
			input.add(KPN(k, "TTTATTT", 6 + 4*i, 6 + 4*i, false)); // 4 + 3*i
			input.add(KPN(k, "AAAGAAA", 6 + 4*i, 6 + 4*i, false)); // 5 + 3*i
			input.add(KPN(k, "AAACAAA", 6 + 4*i, 6 + 4*i, false)); // 6 + 3*i
			KmerPathNode.addEdge(input.get(2 + 3*i), input.get(5 + 3*i));
			KmerPathNode.addEdge(input.get(2 + 3*i), input.get(6 + 3*i));
			KmerPathNode.addEdge(input.get(3 + 3*i), input.get(5 + 3*i));
			KmerPathNode.addEdge(input.get(3 + 3*i), input.get(6 + 3*i));
			KmerPathNode.addEdge(input.get(1 + 3*i), input.get(4 + 3*i));
		}
		// make the leaf longer than all the paths
		// no collapse possible but exhaustive execution required
		input.add(KPN(k, "AAACAAA", 6 + 4*loops, 6 + 4*loops, false));
		KmerPathNode.addEdge(input.get(1 + 3*loops), input.get(4 + 3*loops));
		List<KmerPathNode> result = go(k, 1000000, 1000000, input);
		assertEquals(5 + 3*loops, result.size());
	}
}
