package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.graph.PathNodeBaseFactory;
import au.edu.wehi.idsv.visualisation.NontrackingSubgraphTracker;


public class SubgraphPathNodeTest extends TestHelper {
	@Test
	public void containsReferenceKmer() {
		// no ref kmers if anchor too short
		DeBruijnReadGraph g = RG(4);
		g.addEvidence(SCE(FWD, withSequence("TGTTAAGAC", Read(0, 10, "4M5S"))));
		assertTrue(g.isReference(g.getKmer(KmerEncodingHelper.picardBaseToEncoded(g.getK(), B("TGTT")))));
		assertFalse(g.isReference(g.getKmer(KmerEncodingHelper.picardBaseToEncoded(g.getK(), B("GTTA")))));
		
		DeBruijnPathGraph<DeBruijnSubgraphNode, PathNode<DeBruijnSubgraphNode>> pg = new DeBruijnPathGraph<DeBruijnSubgraphNode, PathNode<DeBruijnSubgraphNode>>(g, new PathNodeBaseFactory<DeBruijnSubgraphNode>(), new NontrackingSubgraphTracker<DeBruijnSubgraphNode, PathNode<DeBruijnSubgraphNode>>());
		assertTrue(pg.isReference(pg.getPaths().iterator().next()));
	}
}
