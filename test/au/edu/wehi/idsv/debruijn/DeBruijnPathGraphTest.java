package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;

import java.util.Arrays;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class DeBruijnPathGraphTest extends TestHelper {
	public class Graph extends DeBruijnGraphBase<DeBruijnNodeBase> {
		public Graph(int k) {
			super(k);
		}
		@Override
		protected DeBruijnNodeBase merge(DeBruijnNodeBase node, DeBruijnNodeBase toAdd) {
			node.add(toAdd);
			return node;
		}
		@Override
		protected DeBruijnNodeBase remove(DeBruijnNodeBase node, DeBruijnNodeBase toRemove) {
			throw new RuntimeException("NYI");
		}
		public Graph add(String sequence) {
			int k = getK();
			assert(sequence.length() >= k);
			byte[] qual = new byte[sequence.length()];
			Arrays.fill(qual, (byte)1);
			for (ReadKmer rk : new ReadKmerIterable(k, B(sequence), qual, false, false)) {
				add(rk.kmer, new DeBruijnNodeBase(rk.weight, new SAMRecord(null)));
			}
			return this;
		}
	}
	public DeBruijnGraphBase<DeBruijnNodeBase> G(int k) {
		return new Graph(k);
	}
	public DeBruijnPathGraph<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>> PG(Graph g, long seed) {
		return new DeBruijnPathGraph<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>>(g, seed, new PathNodeBaseFactory());
	}
	public DeBruijnPathGraph<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>> PG(Graph g) {
		return PG(g, g.getAllKmers().iterator().next());
	}
	@Test
	public void constructor_should_compress_graph_at_construction() {
		DeBruijnPathGraph<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>> pg = PG(new Graph(4)
				.add("GTACCTA"));
		assertEquals(1, pg.getPaths().size());
	}
	@Test
	public void constructor_should_not_compress_forks() {
		DeBruijnPathGraph<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>> pg = PG(new Graph(4)
			.add("GTACCTA")
			.add("GTACCTC"));
		assertEquals(3, pg.getPaths().size());
	}
	@Test
	public void removeSelfIntersectingPaths_should_remove_edges_from_node_to_self() {
		DeBruijnPathGraph<DeBruijnNodeBase, PathNode<DeBruijnNodeBase>> pg = PG(new Graph(4)
			.add("GTACGTA"));
		PathNode<DeBruijnNodeBase> n = pg.getPaths().iterator().next();
		assertEquals(1, pg.nextPath(n).size());
		assertEquals(n, pg.nextPath(n).get(0));
		pg.removeSelfIntersectingPaths();
		assertEquals(0, pg.nextPath(n).size());
	}
}
