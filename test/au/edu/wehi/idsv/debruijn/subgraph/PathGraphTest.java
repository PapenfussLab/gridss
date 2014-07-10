package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;


public class PathGraphTest extends TestHelper {
	@Test
	public void splitOutReferencePaths_should_break_paths_at_reference_non_reference_transition() {
		DeBruijnReadGraph g = G(4, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence("GGTTAACC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(SCE(FWD, withSequence("TTGGT", Read(0, 10, "4M1S"))));
		//     TTAACCGGCCAATT
		//   GGTTAACC
		// TTGGT
		// ^   ^^^^ <- starts of reference kmers
		PathGraph pg = new PathGraph(g, KmerEncodingHelper.picardBaseToEncoded(4, B("TTAA")));
		assertEquals("precondition", 1, pg.getPaths().size());
		pg.splitOutReferencePaths();
		assertEquals(4, pg.getPaths().size());
	}
}
