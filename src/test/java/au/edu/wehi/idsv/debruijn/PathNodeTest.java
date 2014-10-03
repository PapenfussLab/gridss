package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

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
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATACAG"), g);
		assertEquals(6, n.getWeight());
	}
	@Test
	public void getMaxKmerWeight_should_be_largest_kmer_weights() {
		BaseGraph g = G(4)
		.add("ATAC", 1)
		.add("TACA", 2)
		.add("ACAG", 3);
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATACAG"), g);
		assertEquals(3, n.getMaxKmerWeight());
	}
	@Test
	public void length_should_kmer_length() {
		BaseGraph g = G(4)
				.add("ATAG", 4)
				.add("TAGA", 5)
				.add("AGAG", 1)
				.add("GAGT", 8)
				.add("AGTT", 9);
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAGAGTT"), g);
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
		PathNode<DeBruijnNodeBase> unsplit = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAGAGTT"), g);
		// kmer length is 5
		assertEquals(2, new PathNode<DeBruijnNodeBase>(unsplit, 0, 2, g).length());
		assertEquals("ATAGAG", g.S(new PathNode<DeBruijnNodeBase>(unsplit, 0, 3, g)));
		assertEquals("TAGAGT", g.S(new PathNode<DeBruijnNodeBase>(unsplit, 1, 3, g)));
		assertEquals("AGAGTT", g.S(new PathNode<DeBruijnNodeBase>(unsplit, 2, 3, g)));
		assertEquals("AGTT", g.S(new PathNode<DeBruijnNodeBase>(unsplit, 4, 1, g)));
		
		assertEquals(10, new PathNode<DeBruijnNodeBase>(unsplit, 0, 3, g).getWeight());
		assertEquals(14, new PathNode<DeBruijnNodeBase>(unsplit, 1, 3, g).getWeight());
		assertEquals(18, new PathNode<DeBruijnNodeBase>(unsplit, 2, 3, g).getWeight());
		
		assertEquals(5, new PathNode<DeBruijnNodeBase>(unsplit, 0, 3, g).getMaxKmerWeight());
		assertEquals(8, new PathNode<DeBruijnNodeBase>(unsplit, 1, 3, g).getMaxKmerWeight());
		assertEquals(9, new PathNode<DeBruijnNodeBase>(unsplit, 2, 3, g).getMaxKmerWeight());
	}
	@Test
	public void should_merge_shorter_path() {
		BaseGraph g = G(4)
		.add("ATAG", 4)
		.add("TAGA", 5)
		.add("AGAG", 1)
		.add("GAGT", 8)
		.add("AGTT", 9);
		PathNode<DeBruijnNodeBase> unsplit = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAGAGTT"), g);
		
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(g, ImmutableList.of(
			new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAG"), g),
			new PathNode<DeBruijnNodeBase>(toKmer(g, "TAGAG"), g),
			new PathNode<DeBruijnNodeBase>(toKmer(g, "GAGTT"), g)));
		
		assertEquals(unsplit.length(), n.length());
		assertEquals(g.S(unsplit), g.S(n));
		assertEquals(unsplit.getWeight(), n.getWeight());
		assertEquals(unsplit.getMaxKmerWeight(), n.getMaxKmerWeight());
	}
	@Test
	public void should_merge_alt_kmers() {
		BaseGraph g = G(4)
		.add("ATAG", 4)
		.add("TAGA", 5)
		.add("AGAG", 1)
		.add("GAGT", 8)
		.add("AGTT", 9);
		PathNode<DeBruijnNodeBase> unsplit = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAGAGTT"), g);
		
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(g, ImmutableList.of(
			new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAG"), g),
			new PathNode<DeBruijnNodeBase>(toKmer(g, "TAGAG"), g),
			new PathNode<DeBruijnNodeBase>(toKmer(g, "GAGTT"), g)));
		
		assertEquals(unsplit.length(), n.length());
		assertEquals(g.S(unsplit), g.S(n));
		assertEquals(unsplit.getWeight(), n.getWeight());
		assertEquals(unsplit.getMaxKmerWeight(), n.getMaxKmerWeight());
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
		
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAGAGTT"), g);
		n.merge(ImmutableList.of(
				new PathNode<DeBruijnNodeBase>(toKmer(g, "TATCTCA"), g),
				new PathNode<DeBruijnNodeBase>(toKmer(g, "TCAA"), g))
				, g);
		
		assertEquals("ATAGAGTT", g.S(n)); // merge doesn't change main path
		assertEquals(12, n.getMaxKmerWeight()); // merge doesn't change main path
		assertEquals(4+5+1+8+9+2+3+11+1+1, n.getWeight()); // merge doesn't change main path
		assertEquals(10, Lists.newArrayList(Iterables.concat(n.getPathAllKmers())).size());
	}
	@Test
	public void getPathAllKmers_should_include_main_path() {
		BaseGraph g = G(4).add("ATAGAGTT");
		PathNode<DeBruijnNodeBase> n = new PathNode<DeBruijnNodeBase>(toKmer(g, "ATAGAGTT"), g);
		assertEquals(5, Lists.newArrayList(Iterables.concat(n.getPathAllKmers())).size());
	}
}
