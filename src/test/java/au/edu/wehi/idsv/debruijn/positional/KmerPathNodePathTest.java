package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;


public class KmerPathNodePathTest extends TestHelper {
	//  0 - 5   6       8
	//   \   \ / \     /
	//    1 - 2 - 3 - 7 - 9 - 10
	//       /             
	//      4
	//
	public static List<KmerPathNode> testgraph() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "AAAA", 1, 10, false));  // 0
		input.add(KPN(k, "AAATC", 2, 20, false)); // 1
		input.add(KPN(k, "ATCT", 4, 15, false));  // 2
		input.add(KPN(k, "TCTAAT", 1, 50, false));  // 3
		input.add(KPN(k, "GGGATC", 0, 15, false));  // 4
		input.add(KPN(k, "AAAGC", 2, 11, false));  // 5
		input.add(KPN(k, "TCTGAT", 5, 30, false));  // 6
		input.add(KPN(k, "AATCC", 40, 60, false));  // 7
		input.add(KPN(k, "TCCTTATAAATAAC", 42, 62, false));  // 8
		input.add(KPN(k, "TCCGG", 1, 100, false));  // 9
		input.add(KPN(k, "CGGA", 50, 50, false));  // 10
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(5));
		KmerPathNode.addEdge(input.get(1), input.get(2));
		KmerPathNode.addEdge(input.get(5), input.get(2));
		KmerPathNode.addEdge(input.get(4), input.get(2));
		KmerPathNode.addEdge(input.get(2), input.get(6));
		KmerPathNode.addEdge(input.get(2), input.get(3));
		KmerPathNode.addEdge(input.get(6), input.get(3));
		KmerPathNode.addEdge(input.get(3), input.get(7));
		KmerPathNode.addEdge(input.get(7), input.get(8));
		KmerPathNode.addEdge(input.get(7), input.get(9));
		KmerPathNode.addEdge(input.get(9), input.get(10));
		return input;
	}
	public static List<KmerPathNode> sortedKPN(List<KmerPathNode> list) {
		ArrayList<KmerPathNode> out = new ArrayList<KmerPathNode>(list);
		Collections.sort(out, KmerNodeUtil.ByFirstStart);
		return out;
	}
	@Test
	public void dfs_should_iterate_in_start_position_order() {
		List<KmerPathNode> input = testgraph();
		KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(3)), false, 1000);
		assertEquals(1, path.currentPath().size());
		assertTrue(path.dfsNextChild());
		assertEquals(2, path.currentPath().size());
		assertEquals(input.get(2), path.headPath());
		path.dfsPop();
		assertTrue(path.dfsNextChild());
		assertEquals(input.get(6), path.headPath());
		path.dfsPop();
		assertFalse(path.dfsNextChild());
	}
	@Test
	public void asSubnodeList_should_set_consistent_intervals() {
		List<KmerPathNode> input = testgraph();
		for(boolean dir : new boolean[] { true, false }) {
			for (int i = 0; i < input.size(); i++) {
				KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(i)), dir, 1000);
				while (path.dfsNextChild());
				List<KmerPathSubnode> snl = new ArrayList<KmerPathSubnode>(path.headNode().asSubnodes());
				for (int j = 1; j < snl.size(); j++) {
					assertEquals(snl.get(j - 1).firstStart() + snl.get(j - 1).length(), snl.get(j).firstStart());
					assertEquals(snl.get(j - 1).firstEnd() + snl.get(j - 1).length(), snl.get(j).firstEnd());
				}
			}
		}
	}
	@Test
	public void dfs_should_restrict_traversal_interval() {
		List<KmerPathNode> input = testgraph();
		KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(7)), true, 1000);
		assertTrue(path.dfsNextChild());
		assertEquals(input.get(9), path.headPath());
		assertEquals(40, path.headNode().startPositionOfAnchorKmer());
		assertEquals(60, path.headNode().endPositionOfAnchorKmer());
		assertEquals(21, path.headNode().asSubnodes().getFirst().width()); // 9
		assertTrue(path.dfsNextChild());
		assertEquals(input.get(10), path.headPath());
		assertEquals(1, path.headNode().asSubnodes().getFirst().width()); // 10
		assertFalse(path.dfsNextChild());
		path.dfsPop();
		assertFalse(path.dfsNextChild());
		path.dfsPop();
		assertTrue(path.dfsNextChild());
		assertEquals(input.get(8), path.headPath());
		assertEquals(21, path.headNode().asSubnodes().getFirst().width()); // 8
		assertFalse(path.dfsNextChild());
		path.dfsPop();
	}
	@Test
	public void should_limit_traversal_to_max_path_length_of_entire_nodes() {
		List<KmerPathNode> input = testgraph();
		KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(0)), true, 2);
		// no children since adding either 2 or 5 would result in a path length of 3
		assertFalse(path.dfsNextChild());
	}
	@Test
	public void should_list_terminal_nodes_and_leaves() {
		int k = 4;
		List<KmerPathNode> input = new ArrayList<KmerPathNode>();
		input.add(KPN(k, "TAAAAA", 1, 10, false));
		input.add(KPN(k, "AAATC", 5, 5, false));
		input.add(KPN(k, "AAAGC", 8, 9, false));
		input.add(KPN(k, "GCTAA", -1, 7, false));
		input.add(KPN(k, "TAAC", 4, 5, false)); // causes terminal non-leaf since path branches
		input.add(KPN(k, "TAAC", 8, 8, false)); // causes terminal non-leaf since path branches
		input.add(KPN(k, "TGCT", -100, 100, false));
		KmerPathNode.addEdge(input.get(0), input.get(1));
		KmerPathNode.addEdge(input.get(0), input.get(2));
		KmerPathNode.addEdge(input.get(3), input.get(0));
		KmerPathNode.addEdge(input.get(3), input.get(4));
		KmerPathNode.addEdge(input.get(3), input.get(5));
		KmerPathNode.addEdge(input.get(6), input.get(3)); // so root is a leaf
		KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(3)), true, 100);
		assertTrue(path.dfsNextChild());
		assertEquals(input.get(0), path.headPath());
		RangeSet<Integer> expected = TreeRangeSet.create();
		expected.add(Range.closed(1, 1));
		expected.add(Range.closed(3, 4));
		expected.add(Range.closed(7, 9));
		assertEquals(expected, path.headNode().terminalRanges());
		
		expected = TreeRangeSet.create();
		expected.add(Range.closed(-1, -1));
		expected.add(Range.closed(1, 1));
		expected.add(Range.closed(5, 5));
		expected.add(Range.closed(7, 7));
		assertEquals(expected, path.headNode().terminalLeafAnchorRanges());
		
		assertEquals(1, path.headNode().firstTerminalLeaf().node().firstStart());
	}
	@Test
	public void greedyTraverse_should_follow_highest_weighted_path() {
		List<KmerPathNode> input = testgraph();
		KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(0)), true, 100);
		path.greedyTraverse(true, true);
		assertEquals(4, path.headNode().asSubnodes().size());
	}
	@Test
	public void greedyTraverse_should_follow_allowable_nodes() {
		List<KmerPathNode> input = testgraph();
		KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(0)), true, 100);
		path.greedyTraverse(true, false);
		assertEquals(1, path.headNode().asSubnodes().size());
	}
	@Test
	public void greedyTraverse_should_exhaust_alternate_paths() {
		List<KmerPathNode> input = testgraph();
		KmerPathNodePath path = new KmerPathNodePath(new KmerPathSubnode(input.get(0)), true, 100);
		path.greedyTraverse(true, true);
		while (path.headNode() != path.rootNode()) {
			assertFalse(path.dfsNextChild());
			path.dfsPop();
		}
	}
}