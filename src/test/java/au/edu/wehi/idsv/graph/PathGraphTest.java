package au.edu.wehi.idsv.graph;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;

import com.google.common.collect.Lists;

public class PathGraphTest extends TestHelper {
	@Test
	public void constructor_should_compress_graph_at_construction() {
		BasePathGraph pg = PG(G(4)
				.add("GTACCTA"));
		assertEquals(1, pg.getPathCount());
	}
	@Test
	public void constructor_should_not_compress_forks() {
		BasePathGraph pg = PG(G(4)
			.add("GTACCTA")
			.add("GTACCTC"));
		assertEquals(3, pg.getPathCount());
	}
	@Test
	public void ByMaxNodeWeightDesc_should_sort_by_weight_of_kmer() {
		BasePathGraph pg = PG(G(4)
		//     TTAT(GTC)-(GTC)AGTA
		//            /
		//       C(GTC)
				.add("GTCAGTA", 1)
				.add("GTCAGT", 1)
				.add("GTCAG", 1)
				.add("GTCA", 1) // = 4
				.add("TTATGTC", 1)
				.add( "TATGTC", 1)
				.add(  "ATGTC", 1)
				.add(   "TGTC", 3) // = 6
				.add("TCCGTC", 1)
				.add("CCGTC", 1)
				.add("CGTC", 3)); // = 5
		List<PathNode<DeBruijnNodeBase>> nodes = Lists.newArrayList(pg.getPaths());
		assertEquals(3, nodes.size()); // precondition
		Collections.sort(nodes, pg.ByMaxNodeWeightDesc);
		assertEquals(pg.get("TGTC") , nodes.get(0));
		assertEquals(pg.get("CGTC") , nodes.get(1));
		assertEquals(pg.get("GTCA") , nodes.get(2));
	}
}
