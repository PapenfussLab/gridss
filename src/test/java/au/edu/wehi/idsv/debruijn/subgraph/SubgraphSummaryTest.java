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
	@Test
	public void should_track_bounds() {
		SubgraphSummary a = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 10, "10M", true)))), 1, new ReadKmer(0, 2, false));
		DeBruijnSubgraphNode n2 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 20, "10M", true)))), 1, new ReadKmer(0, 2, false));
		a.addNode(n1);
		assertEquals(19, a.getMinAnchor());
		assertEquals(19, a.getMaxAnchor());		
		a.addNode(n2);
		assertEquals(19, a.getMinAnchor());
		assertEquals(29, a.getMaxAnchor());
	}
	@Test
	public void should_merge_bounds() {
		SubgraphSummary a = new SubgraphSummary(0);
		SubgraphSummary b = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 10, "10M", true)))), 1, new ReadKmer(0, 2, false));
		DeBruijnSubgraphNode n2 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 20, "10M", true)))), 1, new ReadKmer(0, 2, false));
		a.addNode(n1);
		b.addNode(n2);
		a.add(b);
		assertEquals(19, a.getMinAnchor());
		assertEquals(29, a.getMaxAnchor());
	}
	@Test
	public void should_track_counts() {
		SubgraphSummary a = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 10, "10M", true)))), 1, new ReadKmer(0, 2, false));
		DeBruijnSubgraphNode n2 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 20, "10M", true)))), 1, new ReadKmer(0, 2, false));
		a.addNode(n1);
		assertEquals(1, a.getKmerCount());		
		a.addNode(n2);
		assertEquals(2, a.getKmerCount());
	}
	@Test
	public void should_merge_counts() {
		SubgraphSummary a = new SubgraphSummary(0);
		SubgraphSummary b = new SubgraphSummary(0);
		DeBruijnSubgraphNode n1 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 10, "10M", true)))), 1, new ReadKmer(0, 2, false));
		DeBruijnSubgraphNode n2 = new DeBruijnSubgraphNode(VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 20, "10M", true)))), 1, new ReadKmer(0, 2, false));
		a.addNode(n1);
		b.addNode(n2);
		a.add(b);
		assertEquals(2, a.getKmerCount());
	}
}
