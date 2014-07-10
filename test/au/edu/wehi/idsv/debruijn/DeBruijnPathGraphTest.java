package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;

import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;


public class DeBruijnPathGraphTest extends TestHelper {
	@Test
	public void constructor_should_compress_graph_at_construction() {
		BasePathGraph pg = PG(G(4)
				.add("GTACCTA"));
		assertEquals(1, pg.getPaths().size());
	}
	@Test
	public void constructor_should_not_compress_forks() {
		BasePathGraph pg = PG(G(4)
			.add("GTACCTA")
			.add("GTACCTC"));
		assertEquals(3, pg.getPaths().size());
	}
	@Test
	public void removeSelfIntersectingPaths_should_remove_edges_from_node_to_self() {
		BasePathGraph pg = PG(G(4)
			.add("GTACGTA"));
		PathNode<DeBruijnNodeBase> n = pg.getPaths().iterator().next();
		assertEquals(1, pg.nextPath(n).size());
		assertEquals(n, pg.nextPath(n).get(0));
		pg.removeSelfIntersectingPaths();
		assertEquals(0, pg.nextPath(n).size());
	}
	@Test
	public void ByMaxKmerWeightDesc_should_sort_by_weight_of_kmer() {
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
				.add("TCCGT", 1)
				.add("CGTC", 3)); // = 5
		List<PathNode<DeBruijnNodeBase>> nodes = Lists.newArrayList(pg.getPaths());
		assertEquals(3, nodes.size()); // precondition
		nodes.sort(pg.ByMaxKmerWeightDesc);
		assertEquals(pg.get("TGTC") , nodes.get(0));
		assertEquals(pg.get("CGTC") , nodes.get(1));
		assertEquals(pg.get("GTCA") , nodes.get(2));
	}
	@Test
	public void ByMaxKmerWeightDesc_should_include_alternate_kmers_in_weight() {
		BasePathGraph pg = PG(G(4)
		//       A(GTC)      
		//            \
		//	TTAT(GTC)-(GTC)AGTA
		//            /
		//       C(GTC)
				.add("AGTC", 2)
				.add("GTCAGTA", 1)
				.add("GTCAGT", 1)
				.add("GTCAG", 1)
				.add("GTCA", 1) // = 4
				.add("TTATGTC", 1)
				.add( "TATGTC", 1)
				.add(  "ATGTC", 1)
				.add(   "TGTC", 3) // = 6
				.add("CGTC", 5)); // = 5 + 2
		pg.mergePaths(ImmutableList.of(pg.get("AGTC")), ImmutableList.of(pg.get("CGTC")));
		List<PathNode<DeBruijnNodeBase>> nodes = Lists.newArrayList(pg.getPaths());
		assertEquals(3, nodes.size()); // precondition
		nodes.sort(pg.ByMaxKmerWeightDesc);
		assertEquals(pg.get("CGTC") , nodes.get(0));
		assertEquals(pg.get("TGTC") , nodes.get(1));
		assertEquals(pg.get("GTCA") , nodes.get(2));
	}
	/**
	 * 
	 */
	@Test
	public void random_sequence_should_not_have_any_repeated_14mers() {
		BasePathGraph pg = PG(G(14).add(S(RANDOM)));
		assertEquals("Other test cases rely on all no branches in this contig", 1, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_should_collapse_resultant_branchless_paths() {
		BasePathGraph pg = PG(G(20)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTGGACGCCAAGTCAGCTGTAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC"));
				//                              ^                ^
		//       *
		//      / \  
		//     *   *        no bubbles
		//      \ /
		//       *
		assertEquals("precondition", 4, pg.getPaths().size());
		pg.collapseSimilarPaths(2, false);
		assertEquals(4, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_should_collapse_only_sufficently_similar_paths() {
		BasePathGraph pg = PG(G(20)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTGGACGCCAAGTCAGCTGTAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC"));
				//                              ^                ^
		//       *
		//      / \  
		//     *   *        
		//      \ /
		//       *
		assertEquals("precondition", 4, pg.getPaths().size());
		pg.collapseSimilarPaths(1, false); // two bases difference = don't collapse
		assertEquals(4, pg.getPaths().size());
		pg.collapseSimilarPaths(2, false); // now we should collapse
		assertEquals(1, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_collapseBubbles_should_not_collapse_non_bubble_branches() {
		BasePathGraph pg = PG(G(16)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGTTAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				//                                   ^^
				.add(                               "TTAAGTCAGCTGAAGCTTTTTTT")
				.add("       AAAAAAAAAAGGTTGTATTCGACGCC")
				);
		//       * - *
		//      / \  
		//     *   *        no bubbles
		//      \ /
		//   * - *
		int before = pg.getPaths().size();
		pg.collapseSimilarPaths(10, false);
		assertEquals(before, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_should_not_collapse_indels() {
		BasePathGraph pg = PG(G(20)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC"));
				//                                              ^ insertion
		assertEquals("precondition", 4, pg.getPaths().size());
		pg.collapseSimilarPaths(10, false);
		assertEquals(4, pg.getPaths().size());
	}
	@Test
	public void mergePaths_should_merge_into_highest_weight_path() {
		BasePathGraph pg = PG(G(4)
		//       A(GTC)      
		//            \
		//	TTAT(GTC)-(GTC)AGTA
		//            /
		//       C(GTC)
				.add("AGTC", 2)
				.add("GTCAGTA", 1)
				.add("GTCAGT", 1)
				.add("GTCAG", 1)
				.add("GTCA", 1) // = 4
				.add("TTATGTC", 1)
				.add( "TATGTC", 1)
				.add(  "ATGTC", 1)
				.add(   "TGTC", 3) // = 6
				.add("CGTC", 5)); // = 5 + 2
		pg.mergePaths(ImmutableList.of(pg.get("AGTC")), ImmutableList.of(pg.get("CGTC")));
		assertNull(pg.get("AGTC"));
		assertEquals(7, pg.get("CGTC").getWeight());
	}
	@Test
	public void shrink_should_remove_self_intersecting_edges() {
		BasePathGraph pg = PG(G(4).add("AAAAATAAAAAA"));
		assertTrue(pg.nextPath(pg.get("AAAA")).contains(pg.get("AAAA")));
		pg.shrink();
		assertFalse(pg.nextPath(pg.get("AAAA")).contains(pg.get("AAAA")));
	}
	@Test
	public void basesDifferent_should_count_bases_on_path() {
		BasePathGraph pg = PG(G(16)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGTTAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				//                                   ^^
				);
		assertEquals(2, pg.basesDifferent(ImmutableList.of(pg.get("TCGACGTTAAGTCAGC")), ImmutableList.of(pg.get("TCGACGCCAAGTCAGC"))));
		String[] s = new String[] {  S(RANDOM).substring(0, 100), S(RANDOM).substring(200, 100) };
		pg = PG(G(16).add(s[0]).add(s[1]));
		int mismatches = 0;
		for (int i = 0; i < 100; i++) if (s[0].charAt(1) != s[1].charAt(1)) mismatches++;
		assertEquals(mismatches, pg.basesDifferent(ImmutableList.of(pg.get(s[0].substring(0, 16))), ImmutableList.of(pg.get(s[1].substring(0, 16)))));
	}
}
