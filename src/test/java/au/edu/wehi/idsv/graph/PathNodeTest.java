package au.edu.wehi.idsv.graph;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;


public class PathNodeTest extends TestHelper {
	@Test
	public void weight_should_be_sum_of_kmer_weights() {
		BaseGraph g = G(4)
		.add("ATAC", 1)
		.add("TACA", 2)
		.add("ACAG", 3);
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATACAG"), g);
		assertEquals(6, n.weight());
	}
	@Test
	public void getMaxKmerWeight_should_be_largest_kmer_weights() {
		BaseGraph g = G(4)
		.add("ATAC", 1)
		.add("TACA", 2)
		.add("ACAG", 3);
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATACAG"), g);
		assertEquals(3, n.maxNodeWeight());
	}
	@Test
	public void length_should_kmer_length() {
		BaseGraph g = G(4)
				.add("ATAG", 4)
				.add("TAGA", 5)
				.add("AGAG", 1)
				.add("GAGT", 8)
				.add("AGTT", 9);
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		assertEquals(5, n.length());
	}
	@Test
	public void should_split_longer_path() {
		BaseGraph g = G(4)
		.add("ATAG", 4)
		.add("TAGA", 5)
		.add("AGAG", 1)
		.add("GAGT", 8)
		.add("AGTT", 9);
		PathNode<DeBruijnNodeBase> unsplit = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		// kmer length is 5
		assertEquals(2, new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 0, 2, g).length());
		assertEquals("ATAGAG", g.S(new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 0, 3, g)));
		assertEquals("TAGAGT", g.S(new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 1, 3, g)));
		assertEquals("AGAGTT", g.S(new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 2, 3, g)));
		assertEquals("AGTT", g.S(new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 4, 1, g)));
		
		assertEquals(10, new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 0, 3, g).weight());
		assertEquals(14, new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 1, 3, g).weight());
		assertEquals(18, new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 2, 3, g).weight());
		
		assertEquals(5, new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 0, 3, g).maxNodeWeight());
		assertEquals(8, new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 1, 3, g).maxNodeWeight());
		assertEquals(9, new PathNode<DeBruijnNodeBase>(ImmutableList.of(unsplit), 2, 3, g).maxNodeWeight());
	}
	@Test
	public void concat_should_append_path() {
		BaseGraph g = G(4)
		.add("ATAG", 4)
		.add("TAGA", 5)
		.add("AGAG", 1)
		.add("GAGT", 8)
		.add("AGTT", 9);
		PathNode<DeBruijnNodeBase> unsplit = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		
		PathNode<DeBruijnNodeBase> n = PathNode.concat(ImmutableList.of(
			new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAG"), g),
			new PathNode<DeBruijnNodeBase>(toNodes(g, "TAGAG"), g),
			new PathNode<DeBruijnNodeBase>(toNodes(g, "GAGTT"), g)),
			g);
		
		assertEquals(unsplit.length(), n.length());
		assertEquals(g.S(unsplit), g.S(n));
		assertEquals(unsplit.weight(), n.weight());
		assertEquals(unsplit.maxNodeWeight(), n.maxNodeWeight());
	}
	@Test
	public void should_merge_alt_kmers() {
		BaseGraph g = G(4)
		.add("ATAG", 4)
		.add("TAGA", 5)
		.add("AGAG", 1)
		.add("GAGT", 8)
		.add("AGTT", 9);
		PathNode<DeBruijnNodeBase> unsplit = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		int weight = unsplit.weight();
		int max = unsplit.maxNodeWeight();				
		
		PathNode<DeBruijnNodeBase> withAlt = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		withAlt.merge(ImmutableList.of(
				new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAG"), g),
				new PathNode<DeBruijnNodeBase>(toNodes(g, "TAGAG"), g),
				new PathNode<DeBruijnNodeBase>(toNodes(g, "GAGTT"), g)),
				g);
		unsplit.merge(ImmutableList.of(withAlt), g);
		
		assertEquals(unsplit.length(), withAlt.length());
		assertEquals(g.S(unsplit), g.S(withAlt));
		assertEquals(3 * weight, unsplit.weight());
		assertEquals(3 * max, unsplit.maxNodeWeight());
	}
	@Test
	public void merge_should_add_merged_kmer_evidence() {
		BaseGraph g = G(4)
		.add("ATAG", 4)
		.add("TAGA", 5)
		.add("AGAG", 1)
		.add("GAGT", 8)
		.add("AGTT", 9)
		
		.add("TATC", 2)
		.add("ATCT", 3)
		.add("TCTC", 11)
		.add("CTCA", 1)
		.add("TCAA", 1);
		
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		n.merge(ImmutableList.of(
				new PathNode<DeBruijnNodeBase>(toNodes(g, "TATCTCA"), g),
				new PathNode<DeBruijnNodeBase>(toNodes(g, "TCAA"), g))
				, g);
		
		assertEquals("ATAGAGTT", g.S(n)); // merge doesn't change main path
		assertEquals(12, n.maxNodeWeight()); // merge doesn't change main path
		assertEquals(4+5+1+8+9+2+3+11+1+1, n.weight()); // merge doesn't change main path
		assertEquals(10, Lists.newArrayList(Iterables.concat(n.getPathAllNodes())).size());
	}
	@Test
	public void merge_should_incorporate_offsets() {
		BaseGraph g = G(4)
		.add("ATAG", 4)
		.add("TAGA", 5)
		.add("AGAG", 1)
		.add("GAGT", 8) // <-- added here
		.add("AGTT", 9) // <-- and here
		
		.add("TATC", 2)
		.add("ATCT", 3) // <-- added
		.add("TCTC", 11) // <-- and this one
		.add("CTCA", 1)
		.add("TCAA", 1);
		// ATAGAGTT
		//   TATCTCAA
		//    ^^
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		n.merge(ImmutableList.of(
				new PathNode<DeBruijnNodeBase>(toNodes(g, "TATCTCA"), g),
				new PathNode<DeBruijnNodeBase>(toNodes(g, "TCAA"), g))
				,1,2,3,
				g);		
		assertEquals("ATAGAGTT", g.S(n));
		assertEquals(9 + 11, n.maxNodeWeight());
		assertEquals(4+5+1+8+9 + 3+11, n.weight());
		for (int i = 0; i < n.length(); i++) {
			assertEquals((i == 3 || i == 4) ? 2 : 1, n.getPathAllNodes().get(i).size());
		}
	}
	@Test
	public void getPathAllNodes_should_include_main_path() {
		BaseGraph g = G(4).add("ATAGAGTT");
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g);
		assertEquals(5, Lists.newArrayList(Iterables.concat(n.getPathAllNodes())).size());
	}
}
