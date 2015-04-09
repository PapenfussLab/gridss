package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.visualisation.NontrackingSubgraphTracker;


public class PathGraphTest extends TestHelper {
	@Test
	public void splitOutReferencePaths_should_break_paths_at_reference_non_reference_transition() {
		DeBruijnReadGraph g = RG(4);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence("GGTTAACC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(SCE(FWD, withSequence("TTGGT", Read(0, 10, "4M1S"))));
		//     TTAACCGGCCAATT
		//   GGTTAACC
		// TTGGT
		// ^   ^^^^ <- starts of reference kmers
		PathGraph pg = new PathGraph(g, KmerEncodingHelper.picardBaseToEncoded(4, B("TTAA")), new NontrackingSubgraphTracker());
		assertEquals("precondition", 1, pg.getPathCount());
		pg.splitOutReferencePaths();
		assertEquals(4, pg.getPathCount());
	}
}
