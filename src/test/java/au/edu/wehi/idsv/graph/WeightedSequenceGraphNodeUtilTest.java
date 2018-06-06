package au.edu.wehi.idsv.graph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;


public class WeightedSequenceGraphNodeUtilTest extends TestHelper {
	@Test
	public void kmerTotalWeight_should_return_path_weight() {
		BaseGraph g = G(4)
			.add("ATAG", 4)
			.add("TAGA", 5)
			.add("AGAG", 1)
			.add("GAGT", 8)
			.add("AGTT", 9);
		assertEquals(4+5+1+8, WeightedSequenceGraphNodeUtil.totalWeight(ImmutableList.of(new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGT"), g))));
		assertEquals(4+5+1+8+9, WeightedSequenceGraphNodeUtil.totalWeight(ImmutableList.of(
				new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGT"), g),
				new PathNode<DeBruijnNodeBase>(toNodes(g, "AGTT"), g))));
	}
	@Test
	public void kmerTotalWeight_should_return_subpath_weight() {
		BaseGraph g = G(4)
			.add("ATAG", 4)
			.add("TAGA", 5)
			.add("AGAG", 1)
			.add("GAGT", 8)
			.add("AGTT", 9)
			.add("GTTC", 10);
		assertEquals(5+1, WeightedSequenceGraphNodeUtil.totalWeight(ImmutableList.of(new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGTT"), g)), 1, 2));
		assertEquals(8+9, WeightedSequenceGraphNodeUtil.totalWeight(ImmutableList.of(
				new PathNode<DeBruijnNodeBase>(toNodes(g, "ATAGAGT"), g),
				new PathNode<DeBruijnNodeBase>(toNodes(g, "AGTTC"), g)),
			3, 2));
	}
	@Test
	public void nodeIterator_should_iterate_over_underlying_primary_nodes() {
		BasePathGraph pg = PG(G(4)
				.add("GTAC")
				.add("TGTA")
				.add("CGTA")
				.add("TACCC")
				.add("TACG")
				);
		List<PathNode<DeBruijnNodeBase>> path = Lists.<PathNode<DeBruijnNodeBase>>newArrayList(pg.get("TGTA"), pg.get("GTAC"), pg.get("TACC"));
		assertEquals(4, Iterators.size(PathNode.nodeIterator(path.iterator())));
		Iterator<DeBruijnNodeBase> it = PathNode.nodeIterator(path.iterator());
		assertEquals(pg.getGraph().getKmer(K("TGTA")), it.next());
		assertEquals(pg.getGraph().getKmer(K("GTAC")), it.next());
		assertEquals(pg.getGraph().getKmer(K("TACC")), it.next());
		assertEquals(pg.getGraph().getKmer(K("ACCC")), it.next());
		assertFalse(it.hasNext());
		
		it = PathNode.nodeIterator(path.iterator(), 1);
		assertEquals(pg.getGraph().getKmer(K("GTAC")), it.next());
		assertEquals(pg.getGraph().getKmer(K("TACC")), it.next());
		assertEquals(pg.getGraph().getKmer(K("ACCC")), it.next());
		assertFalse(it.hasNext());
		
		it = PathNode.nodeIterator(path.iterator(), 2);
		assertEquals(pg.getGraph().getKmer(K("TACC")), it.next());
		assertEquals(pg.getGraph().getKmer(K("ACCC")), it.next());
		assertFalse(it.hasNext());
		
		it = PathNode.nodeIterator(path.iterator(), 3);
		assertEquals(pg.getGraph().getKmer(K("ACCC")), it.next());
		
		it = PathNode.nodeIterator(path.iterator(), 4);
		assertFalse(it.hasNext());
	}
}
