package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class SubgraphSummaryTest extends TestHelper {
	@Test
	public void getRoot_should_not_stack_overflow() {
		SubgraphSummary first = new SubgraphSummary(0);
		SubgraphSummary last = null;
		for (int i = 1; i < 50000; i++) {
			last = new SubgraphSummary(i);
			assertTrue(first.add(last));
		}
		assertEquals(first.getRoot(), last.getRoot()); 
	}
	@Test
	public void merge_should_merge_subgraphs() {
		SubgraphSummary first = new SubgraphSummary(0);
		SubgraphSummary last = new SubgraphSummary(1);
		assertTrue(first.add(last));
		assertEquals(first, first.getRoot()); 
		assertEquals(first, last.getRoot());
	}
}
