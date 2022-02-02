package au.edu.wehi.idsv.graph;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import com.google.common.collect.ImmutableList;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


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
}
