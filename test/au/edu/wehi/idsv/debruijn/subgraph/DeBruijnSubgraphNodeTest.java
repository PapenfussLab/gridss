package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;


public class DeBruijnSubgraphNodeTest extends TestHelper {
	@Test
	public void should_set_reference_bounds_based_on_evidence() {
		VariantEvidence e = VariantEvidence.createSoftClipEvidence(FWD, 4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		for (int i = 0; i < 4; i++) {
			assertEquals(e.getInferredReferencePosition(i), (int)new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2)).getMinReferencePosition());
			assertEquals(e.getInferredReferencePosition(i), (int)new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2)).getMaxReferencePosition());
			assertTrue(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2)).isReference());
		}
		for (int i = 4; i < 11; i++) {
			assertNull(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2)).getMinReferencePosition());
			assertNull(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2)).getMaxReferencePosition());
			assertFalse(new DeBruijnSubgraphNode(e, i, new ReadKmer(0, 2)).isReference());
		}
	}
	@Test
	public void should_merge_reference_position_bounds() {
		VariantEvidence e1 = VariantEvidence.createSoftClipEvidence(FWD, 4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		VariantEvidence e2 = VariantEvidence.createSoftClipEvidence(FWD, 4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 11, "7M7S"))));
		for (int i = 0; i < 4; i++) {
			DeBruijnSubgraphNode n1 = new DeBruijnSubgraphNode(e1, i, new ReadKmer(0, 2));
			DeBruijnSubgraphNode n2 = new DeBruijnSubgraphNode(e2, i, new ReadKmer(0, 2));
			n1.add(n2);
			assertEquals(e1.getInferredReferencePosition(i), (int)n1.getMinReferencePosition());
			assertEquals(e2.getInferredReferencePosition(i), (int)n1.getMaxReferencePosition());
		}
	}
}
