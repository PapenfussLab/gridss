package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.visualisation.StaticDeBruijnPathGraphGexfExporter;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSortedSet;
import org.junit.Test;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import static org.junit.Assert.assertEquals;


public class DeBruijnPathGraphTest extends TestHelper {
	@Test
	public void random_sequence_should_not_have_any_repeated_14mers() {
		BasePathGraph pg = PG(G(14).add(S(RANDOM)));
		assertEquals("Other test cases rely on all no branches in this contig", 1, pg.getPathCount());
	}
	// Collapse Leaf test cases:
	//      *   Leaf  *
	//     / \  here / \
	//  - *   *  -  *   * -
	//     \ /       \ /
	//      A         B
	public BasePathGraph leafpathTestPathGraph(String leaf) {
		return leafpathTestPathGraph(leaf, 1);
	}
	public BasePathGraph leafpathTestPathGraph(String leaf, int leafWeight) {
		return PG(G(3)
				.add("AAACGGATCT", 10) // main path
				.add("AAGG") // A
				.add("ATTC") // B
				.add(leaf, leafWeight));
	}
	private void debugDumpToCDevDebugGraphGexf(BasePathGraph pg, int k) {
		StaticDeBruijnPathGraphGexfExporter<DeBruijnNodeBase, DeBruijnPathNode<DeBruijnNodeBase>> tmp =
				new StaticDeBruijnPathGraphGexfExporter<DeBruijnNodeBase, DeBruijnPathNode<DeBruijnNodeBase>>();
		tmp.snapshot(pg);
		tmp.saveTo(new File("C:/dev/debugGraph.gexf"));
	}
	@Test
	public void shrink_should_remove_self_intersecting_path_edges() {
		BasePathGraph pg = PG(G(4)
			.add("GTAC"));
		// lets add a self-intersecting edge
		pg.addEdge(pg.get("GTAC"), pg.get("GTAC"));
		assertEquals(1, pg.next(pg.get("GTAC")).size());
		pg.shrink();
		assertEquals(0, pg.next(pg.get("GTAC")).size());
	}
	
	@Test
	public void seed_should_restrict_path_graph_to_reachable_subgraph() {
		String[] s = new String[] { S(RANDOM).substring(0, 100), S(RANDOM).substring(200, 300) };
		assertEquals(2, PG(G(20).add(s[0]).add(s[1])).getPathCount());
		assertEquals(1, PG(G(20).add(s[0]).add(s[1]), s[0].substring(0, 20)).getPathCount());
	}
	@Test
	public void splitPathToEnsureBreaksAt_should_split_at_given_offsets() {
		BasePathGraph pg = PG(G(4).add("TTAAGGCCGGTTGGAACC"));
		//                                 123456789012345
		assertEquals("precondition", 1, pg.getPathCount());
		assertEquals("precondition", 15, pg.get("TTAA").length());
		List<DeBruijnPathNode<DeBruijnNodeBase>> result = pg.splitPathToEnsureBreaksAt(ImmutableList.of(pg.get("TTAA")), ImmutableSortedSet.of(0, 2, 5, 15));
		assertEquals(3, result.size());
		assertEquals("TTAAG", pg.S(result.get(0)));
		assertEquals("AAGGCC", pg.S(result.get(1)));
		assertEquals("GCCGGTTGGAACC", pg.S(result.get(2)));
		result = pg.splitPathToEnsureBreaksAt(result, ImmutableSortedSet.of(1,4,10,11,12,13));
		// path bounds at: 0,1,2,4,5,10,11,12,13,15
		// length =         1 1 2 1 5  1  1  1  2  
		assertEquals(9, result.size());
		assertEquals(1, result.get(0).length());;
		assertEquals(1, result.get(1).length());
		assertEquals(2, result.get(2).length());
		assertEquals(1, result.get(3).length());
		assertEquals(5, result.get(4).length());
		assertEquals(1, result.get(5).length());
		assertEquals(1, result.get(6).length());
		assertEquals(1, result.get(7).length());
		assertEquals(2, result.get(8).length());
		assertEquals("GAACC", pg.S(result.get(8)));
	}
	@Test
	public void greedyTraverse_should_follow_first_comparator_path() {
		BasePathGraph pg = PG(G(4)
				.add("GTAC")
				.add("TACT", 1)
				.add("TACA", 2)
				.add("TACG", 3)
				.add("TACC", 4)
				.add("AGTA", 1)
				.add("CGTA", 2)
				.add("GGTA", 3)
				.add("TGTA", 4)
				);
		LinkedList<DeBruijnPathNode<DeBruijnNodeBase>> path = pg.greedyTraverse(pg.get("GTAC"), pg.ByPathTotalWeightDesc, pg.ByPathTotalWeightDesc, null);
		assertEquals("TGTACC", pg.S(path));
	}
	@Test
	public void greedyTraverse_should_not_follow_excluded_path() {
		BasePathGraph pg = PG(G(4)
				.add("GTAC")
				.add("TACT", 1)
				.add("TACA", 2)
				.add("TACG", 3)
				.add("TACC", 4)
				);
		LinkedList<DeBruijnPathNode<DeBruijnNodeBase>> path = pg.greedyTraverse(pg.get("GTAC"), pg.ByPathTotalWeightDesc, pg.ByPathTotalWeightDesc, ImmutableSet.of(pg.get("TACC")));
		assertEquals("GTACG", pg.S(path));
	}
}
