package au.edu.wehi.idsv.graph;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNode;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertEquals;

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
		List<DeBruijnPathNode<DeBruijnNodeBase>> nodes = Lists.newArrayList(pg.getPaths());
		assertEquals(3, nodes.size()); // precondition
		Collections.sort(nodes, pg.ByMaxNodeWeightDesc);
		assertEquals(pg.get("TGTC") , nodes.get(0));
		assertEquals(pg.get("CGTC") , nodes.get(1));
		assertEquals(pg.get("GTCA") , nodes.get(2));
	}
	@Test
	public void should_reduce_simple_paths() {
		BasePathGraph pg = PG(G(4)
				.add("AAAATTT"));
		List<DeBruijnPathNode<DeBruijnNodeBase>> nodes = Lists.newArrayList(pg.getPaths());
		assertEquals(1, nodes.size());
	}
	@Test
	public void should_not_reduce_forks() {
		BasePathGraph pg = PG(G(4)
				.add("TCCGAAT")
				.add("TCCGAAC")
				.add("TCCGAAG"));
		List<DeBruijnPathNode<DeBruijnNodeBase>> nodes = Lists.newArrayList(pg.getPaths());
		assertEquals(4, nodes.size());
	}
}
