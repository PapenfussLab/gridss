package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.VariantEvidence;


public class DeBruijnSubgraphNodeTest extends TestHelper {
	private static VariantEvidence fsc = new VariantEvidence(4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))), getContext().getLinear());
	@Test
	public void should_track_bounds() {
		DeBruijnSubgraphNode r0 = new DeBruijnSubgraphNode(fsc, 0);
		DeBruijnSubgraphNode r1 = new DeBruijnSubgraphNode(fsc, 0);
		DeBruijnSubgraphNode r2 = new DeBruijnSubgraphNode(fsc, 0);
		assertEquals(r1.getExpectedPosition(), r1.getMinLinearPosition(), 0);
		assertEquals(r1.getExpectedPosition(), r1.getMaxLinearPosition(), 0);
		r1.add(r0);
		assertEquals(r0.getExpectedPosition(), r1.getMinLinearPosition(), 0);
		assertEquals(r1.getExpectedPosition(), r1.getMaxLinearPosition(), 0);
		r1.add(r2);
		assertEquals(r0.getExpectedPosition(), r1.getMinLinearPosition(), 0);
		assertEquals(r2.getExpectedPosition(), r1.getMaxLinearPosition(), 0);
	}
}
