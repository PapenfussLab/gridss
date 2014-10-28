package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.StaticDeBruijnPathGraphGexfExporter;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Lists;


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
				.add("CCGTC", 1)
				.add("CGTC", 3)); // = 5
		List<PathNode<DeBruijnNodeBase>> nodes = Lists.newArrayList(pg.getPaths());
		assertEquals(3, nodes.size()); // precondition
		Collections.sort(nodes, pg.ByMaxKmerWeightDesc);
		assertEquals(pg.get("TGTC") , nodes.get(0));
		assertEquals(pg.get("CGTC") , nodes.get(1));
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
	public void collapseSimilarPaths_should_collapse_resultant_branchless_paths() throws AlgorithmRuntimeSafetyLimitExceededException {
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
		assertEquals(1, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_should_collapse_only_sufficently_similar_paths() throws AlgorithmRuntimeSafetyLimitExceededException {
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
	public void collapseSimilarPaths_should_merge_dissimilar_nodes() throws AlgorithmRuntimeSafetyLimitExceededException {
		BasePathGraph pg = PG(G(20)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTGGACGCCATGTCAGCTGTAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC"));
				//                              ^       ^        ^
		pg.splitPathToEnsureBreaksAt(ImmutableList.of(pg.get("GGACGCCATGTCAGCTGTAG")), ImmutableSortedSet.of(1,2,3,5,7));
		pg.splitPathToEnsureBreaksAt(ImmutableList.of(pg.get("CGACGCCAAGTCAGCTGAAG")), ImmutableSortedSet.of(4,6,8));
		pg.collapseSimilarPaths(3, false); // now we should collapse
		assertEquals(1, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_collapseBubbles_should_not_collapse_non_bubble_branches() throws AlgorithmRuntimeSafetyLimitExceededException {
		BasePathGraph pg = PG(G(4)
				.add("TTAAGGCCGTGTTGGAACC")
				.add("TTAAGGCCGAGTTGGAACC")
				.add(       "TGAG")
				.add(        "GTGA")
				//             ^
				);
		//            5  - 2
		//		    /  \     \
		//		  1      4    0    = no bubbles
		//		    \     \  /
		//		     6  -  3
		assertEquals(7, pg.getPaths().size());
		pg.collapseSimilarPaths(10, true);
		assertEquals(7, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_should_not_collapse_indels() throws AlgorithmRuntimeSafetyLimitExceededException {
		BasePathGraph pg = PG(G(20)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC"));
				//                                              ^ insertion
		assertEquals("precondition", 4, pg.getPaths().size());
		pg.collapseSimilarPaths(10, false);
		assertEquals(4, pg.getPaths().size());
	}
	@Test
	public void collapseSimilarPaths_should_collapse_branches() throws AlgorithmRuntimeSafetyLimitExceededException {
		BasePathGraph pg = PG(G(20)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACCGATC"));
				//                                                                                        ^^
		assertEquals("precondition", 3, pg.getPaths().size());
		pg.collapseSimilarPaths(1, false);
		assertEquals(3, pg.getPaths().size());
		pg.collapseSimilarPaths(2, false);
		assertEquals(1, pg.getPaths().size());
	}
	@Test
	public void collapseLeaves_should_collapse_branches_of_different_lengths() throws AlgorithmRuntimeSafetyLimitExceededException {
		BasePathGraph pg = PG(G(20)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACCGATCGATCA", 2));
				//                                                                                        ^^
		assertEquals("precondition", 3, pg.getPaths().size());
		pg.collapseSimilarPaths(2, false);
		assertEquals(3, pg.getPaths().size());
		pg.collapseLeaves(2);
		assertEquals(1, pg.getPaths().size());
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
	@Test
	public void collapseLeaves_should_collapse_only_leaves() throws AlgorithmRuntimeSafetyLimitExceededException {
		BasePathGraph pg = leafpathTestPathGraph("AAA");
		assertEquals(0, pg.collapseLeaves(3));
	}
	@Test
	public void collapseLeaves_should_collapse_to_leaf_to_path() throws AlgorithmRuntimeSafetyLimitExceededException {
		// forward
		BasePathGraph pg = leafpathTestPathGraph("GATGT");//"GGGTCA");
		assertEquals("test case precondition", 8, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(0));
		assertEquals(1, pg.collapseLeaves(1));
		assertEquals(7, pg.paths.size());
		// backward
		pg = leafpathTestPathGraph("GTGAT");
		assertEquals("test case precondition", 9, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(1));
		assertEquals(1, pg.collapseLeaves(2));
		assertEquals(7, pg.paths.size());
	}
	private void debugDumpToCDevDebugGraphGexf(BasePathGraph pg, int k) {
		StaticDeBruijnPathGraphGexfExporter<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>> tmp =
				new StaticDeBruijnPathGraphGexfExporter<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>>(k);
		tmp.snapshot(pg);
		tmp.saveTo(new File("C:/dev/debugGraph.gexf"));
	}
	@Test
	public void collapseLeaves_should_not_collapse_into_lower_weight_path() throws AlgorithmRuntimeSafetyLimitExceededException {
		// forward
		BasePathGraph pg = leafpathTestPathGraph("GATGT", 10);
		debugDumpToCDevDebugGraphGexf(pg, 3);
		assertEquals(0, pg.collapseLeaves(3));
		// backward
		pg = leafpathTestPathGraph("GTGAT", 10);
		assertEquals(0, pg.collapseLeaves(3));
	}
	@Test
	public void collapseLeaves_should_collapse_smaller_to_leaf_to_leaf() throws AlgorithmRuntimeSafetyLimitExceededException {
		// forward
		BasePathGraph pg = PG(G(4)
				.add("AAAAT",2)
				.add("AAAAC"));
		assertEquals("test case precondition", 3, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(0));
		assertEquals(1, pg.collapseLeaves(1));
		assertEquals(1, pg.paths.size());
		// Should merge the path then collapse the two remaining nodes into one
		assertEquals("AAAAT", pg.paths.iterator().next().pathKmerString(4));
		// backward
		pg = PG(G(4)
				.add("TAAAA",2)
				.add("CAAAA"));
		assertEquals("test case precondition", 3, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(0));
		assertEquals(1, pg.collapseLeaves(1));
		assertEquals(1, pg.paths.size());
		assertEquals("TAAAA", pg.paths.iterator().next().pathKmerString(4));
	}
	@Ignore("See github issue 12: de Bruijn subgraph contig construction should call highest weighted base")
	@Test
	public void collapseLeaves_leaf_leaf_collapse_should_collapse_based_to_kmer_weight_of_shared_bases() throws AlgorithmRuntimeSafetyLimitExceededException {
		BasePathGraph pg = PG(G(4)
				.add("AAAAT",2)
				.add("AAAAC")
				.add("AACCC", 3));
		assertEquals("test case precondition", 3, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(0));
		assertEquals(1, pg.collapseLeaves(1));
		assertEquals(1, pg.paths.size());
		// take the AAAAT path then extend with the additional bases from the other path
		assertEquals("AAAATCC", pg.paths.iterator().next().pathKmerString(4));
		// backward
		pg = PG(G(4)
				.add("TAAAA",2)
				.add("CAAAA")
				.add("CCCAA", 3));
		assertEquals("test case precondition", 3, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(0));
		assertEquals(1, pg.collapseLeaves(1));
		assertEquals(1, pg.paths.size());
		assertEquals("CCTGAAAA", pg.paths.iterator().next().pathKmerString(4));
	}
	/**
	 * See github issue 12: de Bruijn subgraph contig construction should call highest weighted base
	 * https://github.com/d-cameron/idsv/issues/12
	 */
	@Test
	public void collapseLeaves_should_not_collapse_longer_leaf()  throws AlgorithmRuntimeSafetyLimitExceededException{
		BasePathGraph pg = PG(G(4)
				.add("AAAAT",2)
				.add("AAAAC")
				.add("AACCC", 3));
		assertEquals("test case precondition", 3, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(8));
		pg = PG(G(4)
				.add("TAAAA",2)
				.add("CAAAA")
				.add("CCCAA", 3));
		assertEquals("test case precondition", 3, pg.paths.size());
		assertEquals(0, pg.collapseLeaves(8));
	}
	@Ignore("See github issue 12: de Bruijn subgraph contig construction should call highest weighted base")
	@Test
	public void collapseLeaves_should_fix_inconsistent_kmers()  throws AlgorithmRuntimeSafetyLimitExceededException{
		// leaf-leaf collapse results in a main kmer path that doesn't
		// make sense.
		// TODO: FIXME: replace the inconsistent kmers with the correct ones.
		BasePathGraph pg = PG(G(4)
				.add("AAAAT",2)
				.add("AAAAC")
				.add("AACC"));
		pg.collapseLeaves(1);
		// if we just merged the kmers then the main path would be inconsistent
		// as it would go
		// AAAA -> AAAT -> AACC
		//              ^     \___ need to change this C to a T
		//              |
		//         bad transition
		PathNode<DeBruijnNodeBase> pn = pg.paths.iterator().next();
		assertEquals("AATC", S(KmerEncodingHelper.encodedToPicardBases(4, pn.getLast())));
		
		// test backwards as well
		pg = PG(G(4)
				.add("TAAAA",2)
				.add("CAAAA")
				.add("CCAA"));
		pg.collapseLeaves(1);
		pn = pg.paths.iterator().next();
		assertEquals("CTAA", S(KmerEncodingHelper.encodedToPicardBases(4, pn.getFirst())));
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
		assertEquals(pg.get("CGTC"), pg.get("AGTC"));
		assertEquals(7, pg.get("CGTC").getWeight());
	}
	@Test
	public void shrink_should_remove_self_intersecting_path_edges() {
		BasePathGraph pg = PG(G(4)
			.add("GTAC"));
		// lets add a self-intersecting edge
		pg.addEdge(pg.get("GTAC"), pg.get("GTAC"));
		assertEquals(1, pg.nextPath(pg.get("GTAC")).size());
		pg.shrink();
		assertEquals(0, pg.nextPath(pg.get("GTAC")).size());
	}
	@Test
	public void basesDifferent_should_count_bases_on_path() {
		BasePathGraph pg = PG(G(16)
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				.add("CATTAATCGCAAGAGCGGGTTGTATTCGACGTTAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATC")
				//                                   ^^
				);
		assertEquals(2, pg.basesDifferent(ImmutableList.of(pg.get("TCGACGTTAAGTCAGC")), ImmutableList.of(pg.get("TCGACGCCAAGTCAGC"))));
		
		
		String[] s = new String[] {  S(RANDOM).substring(0, 100), S(RANDOM).substring(200, 300) };
		pg = PG(G(20)
			.add(s[0])
			.add(s[1]));
		int mismatches = 0;
		for (int i = 0; i < 100; i++) if (s[0].charAt(i) != s[1].charAt(i)) mismatches++;
		assertEquals(mismatches, pg.basesDifferent(
				ImmutableList.of(pg.get(s[0].substring(0, 20))),
				ImmutableList.of(pg.get(s[1].substring(0, 20)))));
	}
	@Test
	public void seed_should_restrict_path_graph_to_reachable_subgraph() {
		String[] s = new String[] { S(RANDOM).substring(0, 100), S(RANDOM).substring(200, 300) };
		assertEquals(2, PG(G(20).add(s[0]).add(s[1])).getPaths().size());
		assertEquals(1, PG(G(20).add(s[0]).add(s[1]), s[0].substring(0, 20)).getPaths().size());
	}
	@Test
	public void splitPathToEnsureBreaksAt_should_split_at_given_offsets() {
		BasePathGraph pg = PG(G(4).add("TTAAGGCCGGTTGGAACC"));
		//                                 123456789012345
		assertEquals("precondition", 1, pg.getPaths().size());
		assertEquals("precondition", 15, pg.get("TTAA").length());
		List<PathNode<DeBruijnNodeBase>> result = pg.splitPathToEnsureBreaksAt(ImmutableList.of(pg.get("TTAA")), ImmutableSortedSet.of(0, 2, 5, 15));
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
		LinkedList<PathNode<DeBruijnNodeBase>> path = pg.greedyTraverse(pg.get("GTAC"), pg.ByPathTotalWeightDesc, pg.ByPathTotalWeightDesc, null);
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
		LinkedList<PathNode<DeBruijnNodeBase>> path = pg.greedyTraverse(pg.get("GTAC"), pg.ByPathTotalWeightDesc, pg.ByPathTotalWeightDesc, ImmutableSet.of(pg.get("TACC")));
		assertEquals("GTACG", pg.S(path));
	}
	@Test
	public void mergePaths_should_not_merge_inconsistent_paths() {
		BasePathGraph pg = PG(G(1)
				.add("AAGATACAGTGCGGTCTTCC"));
		assertFalse(pg.mergePaths(ImmutableList.of(pg.get("A"), pg.get("B")), ImmutableList.of(pg.get("B"), pg.get("A"))));
	}
	@Test
	public void mergePaths_should_compare_path_position_of_nodes_when_checking_for_inconsistent_paths() {
		//   1 - 1 - B
		//  /
		// A             paths are fine to merge since B is at the same position on both paths
		//  \
		//   2 - B
		BasePathGraph pg = PG(G(4)
				.add("AAAA") // A
				.add("TTTT") // 1
				.add("CCCC") // 1
				.add("GGGGT") // 2
				.add("GTAC") // B
				);
		assertTrue(pg.mergePaths(
				ImmutableList.of(pg.get("AAAA"), pg.get("TTTT"), pg.get("CCCC"), pg.get("GTAC")),
				ImmutableList.of(pg.get("AAAA"), pg.get("GGGG"), pg.get("GTAC"))));
	}
}
