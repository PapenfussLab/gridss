package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;


public class DeBruijnSubgraphNodeTest extends TestHelper {
	@Test
	public void should_set_reference_bounds_based_on_evidence() {
		VariantEvidence e = VariantEvidence.createSoftClipEvidence(FWD, 4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		for (int i = 0; i < 4; i++) {
			assertEquals(e.getExpectedReferencePosition(i), (int)new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2, false)).getMinReferencePosition());
			assertEquals(e.getExpectedReferencePosition(i), (int)new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2, false)).getMaxReferencePosition());
			assertTrue(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2, false)).isReference());
		}
		for (int i = 4; i < 11; i++) {
			assertNull(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2, false)).getMinReferencePosition());
			assertNull(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2, false)).getMaxReferencePosition());
			assertFalse(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2, false)).isReference());
		}
	}
	@Test
	public void should_merge_reference_position_bounds() {
		VariantEvidence e1 = VariantEvidence.createSoftClipEvidence(FWD, 4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		VariantEvidence e2 = VariantEvidence.createSoftClipEvidence(FWD, 4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 11, "7M7S"))));
		for (int i = 0; i < 4; i++) {
			DeBruijnSubgraphNode n1 = new DeBruijnSubgraphNode(e1, i, new ReadKmer(0, 2, false));
			assertFalse(n1.isMateAnchored());
			assertTrue(n1.isReference());
			assertEquals(e1.getExpectedReferencePosition(i), (int)n1.getMinReferencePosition());
			assertEquals(e1.getExpectedReferencePosition(i), (int)n1.getMaxReferencePosition());
			assertNull(n1.getMaxMatePosition());
			assertNull(n1.getMinMatePosition());
			DeBruijnSubgraphNode n2 = new DeBruijnSubgraphNode(e2, i, new ReadKmer(0, 2, false));
			assertFalse(n2.isMateAnchored());
			assertTrue(n2.isReference());
			assertEquals(e2.getExpectedReferencePosition(i), (int)n2.getMinReferencePosition());
			assertEquals(e2.getExpectedReferencePosition(i), (int)n2.getMaxReferencePosition());
			assertNull(n2.getMaxMatePosition());
			assertNull(n2.getMinMatePosition());
			n1.add(n2);
			n1.add(n2);
			assertEquals(e1.getExpectedReferencePosition(i), (int)n1.getMinReferencePosition());
			assertEquals(e2.getExpectedReferencePosition(i), (int)n1.getMaxReferencePosition());
			assertEquals(e2.getExpectedReferencePosition(i), (int)n1.getBestReferencePosition());
		}
	}
	@Test
	public void should_merge_mate_bounds() {
		VariantEvidence e1 = VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 10, "10M", true))));
		VariantEvidence e2 = VariantEvidence.createRemoteReadEvidence(FWD,  4,  NRRP(withSequence("TTAACCGGCC", OEA(0, 20, "10M", true))));
		for (int i = 0; i < 4; i++) {
			DeBruijnSubgraphNode n1 = new DeBruijnSubgraphNode(e1, i, new ReadKmer(0, 2, false));
			assertEquals(19, (int)n1.getMaxMatePosition());
			assertEquals(19, (int)n1.getMinMatePosition());
			assertNull(n1.getMinReferencePosition());
			assertNull(n1.getMaxReferencePosition());
			assertTrue(n1.isMateAnchored());
			assertFalse(n1.isReference());
			DeBruijnSubgraphNode n2 = new DeBruijnSubgraphNode(e2, i, new ReadKmer(0, 2, false));
			assertEquals(29, (int)n2.getMaxMatePosition());
			assertEquals(29, (int)n2.getMinMatePosition());
			assertNull(n2.getMinReferencePosition());
			assertNull(n2.getMaxReferencePosition());
			assertTrue(n2.isMateAnchored());
			assertFalse(n2.isReference());
			n1.add(n2);
			assertEquals(19, (int)n1.getMinMatePosition());
			assertEquals(29, (int)n1.getMaxMatePosition());
			assertNull(n1.getMinReferencePosition());
			assertNull(n1.getMaxReferencePosition());
			assertTrue(n1.isMateAnchored());
			assertFalse(n1.isReference());
		}
	}
}
