package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class TraversalNodeTest extends TestHelper {
	@Test
	public void contains_should_check_overlapping_nodes_only() {
		KmerPathNode n1 = KPN(4, "TTTT", 1, 10, false);
		KmerPathNode n2 = KPN(4, "CCCC", 7, 15, false);
		KmerPathNode n3 = KPN(4, "AAAA", 11, 20, false);
		TraversalNode tn1 = new TraversalNode(new KmerPathSubnode(n1, 1, 10), 0);
		TraversalNode tn2 = new TraversalNode(tn1, new KmerPathSubnode(n2, 11, 12));
		TraversalNode tn3 = new TraversalNode(tn2, new KmerPathSubnode(n3, 12, 12));
		assertTrue(tn1.traversingWouldCauseSelfIntersection(n1));
		assertFalse(tn1.traversingWouldCauseSelfIntersection(n2));
		assertFalse(tn1.traversingWouldCauseSelfIntersection(n3));
		assertTrue(tn2.traversingWouldCauseSelfIntersection(n1));
		assertTrue(tn2.traversingWouldCauseSelfIntersection(n2));
		assertFalse(tn2.traversingWouldCauseSelfIntersection(n3));
		assertTrue(tn3.traversingWouldCauseSelfIntersection(n2));
		assertTrue(tn3.traversingWouldCauseSelfIntersection(n3));
		
		// whilst n1 is in fact on path, it is not possible for
		// n1 to cause self-intersection because adding it to tn3
		// is not possible.
		assertFalse(tn3.traversingWouldCauseSelfIntersection(n1)); // result undefined since traversal is invalid
	}
}
