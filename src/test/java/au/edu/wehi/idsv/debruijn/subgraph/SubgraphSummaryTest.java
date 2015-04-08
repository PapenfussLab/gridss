package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;


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
	private DeBruijnSubgraphNode node(int pos) {
		VariantEvidence e = new VariantEvidence(4, SCE(FWD, Read(0, 1, "500M500S")), getContext().getLinear());
		DeBruijnSubgraphNode n = new DeBruijnSubgraphNode(e, pos - 1,  new ReadKmer(0, 2, false));
		return n;
	}
	@Test
	public void should_track_bounds() {
		SubgraphSummary a = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = node(19);
		DeBruijnSubgraphNode n2 = node(29);
		a.addNode(n1);
		assertEquals(LCCB + 19, a.getMinLinearPosition());
		assertEquals(LCCB + 19, a.getMaxLinearPosition());	
		a.addNode(n2);
		assertEquals(LCCB + 19, a.getMinLinearPosition());
		assertEquals(LCCB + 29, a.getMaxLinearPosition());
	}
	@Test
	public void should_merge_bounds() {
		SubgraphSummary a = new SubgraphSummary(0);
		SubgraphSummary b = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = node(19);
		DeBruijnSubgraphNode n2 = node(29);
		a.addNode(n1);
		b.addNode(n2);
		a.add(b);
		assertEquals(LCCB + 19, a.getMinLinearPosition());
		assertEquals(LCCB + 29, a.getMaxLinearPosition());
	}
	@Test
	public void should_track_counts() {
		SubgraphSummary a = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = node(19);
		DeBruijnSubgraphNode n2 = node(29);
		a.addNode(n1);
		assertEquals(1, a.getKmerCount());		
		a.addNode(n2);
		assertEquals(2, a.getKmerCount());
	}
	@Test
	public void should_merge_counts() {
		SubgraphSummary a = new SubgraphSummary(0);
		SubgraphSummary b = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = node(19);
		DeBruijnSubgraphNode n2 = node(29);
		a.addNode(n1);
		b.addNode(n2);
		a.add(b);
		assertEquals(2, a.getKmerCount());
	}
}
