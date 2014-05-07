package au.edu.wehi.socrates.graph;

import java.util.List;

import org.junit.Test;

import static org.junit.Assert.*;

import com.google.common.collect.Lists;

public class TrapezoidGraphTest {
	private TrapezoidGraphNode N(long start1, long end1, long start2, long end2) {
		return new TrapezoidGraphNode(start1, end1, start2, end2, 1);
	}
	private TrapezoidGraphNode N(long start1, long end1, long start2, long end2, double weight) {
		return new TrapezoidGraphNode(start1, end1, start2, end2, weight);
	}
	TrapezoidGraph graph; 
	private TrapezoidGraphNode[] getCliques(TrapezoidGraphNode[] nodes) {
		graph = new TrapezoidGraph();
		for (TrapezoidGraphNode n : nodes) {
			graph.add(n);
		}
		List<TrapezoidGraphNode> result = Lists.newArrayList(graph.getAllMaximalCliques());
		return result.toArray(new TrapezoidGraphNode[0]);
	}
	private void assertContains(TrapezoidGraphNode[] expected, TrapezoidGraphNode test) {
		for (TrapezoidGraphNode n : expected) {
			if (test.start1 == n.start1
					&& test.end1 == n.end1
					&& test.start2 == n.start2
					&& test.end2 == n.end2
					&& test.weight == n.weight) {
				return;
			}
		}
		fail("Unexpected maximal clique call: " + test.toString());
	}
	private TrapezoidGraphNode[] flipXY(TrapezoidGraphNode[] a) {
		TrapezoidGraphNode[] r = new TrapezoidGraphNode[a.length];
		for (int i = 0; i < a.length; i++) {
			r[i] = new TrapezoidGraphNode(a[i].start2, a[i].end2, a[i].start1, a[i].end1, a[i].weight);
		}
		return r;
	}
	private void assertMatches(TrapezoidGraphNode[] expected, TrapezoidGraphNode[] actual) {
		assertEquals("Number of cliques called does not match expected number", expected.length, actual.length);
		for (TrapezoidGraphNode n : actual) {
			assertContains(expected, n);
		}
	}
	private void go(TrapezoidGraphNode[] input, TrapezoidGraphNode... expected) {		
		assertMatches(expected, getCliques(input));
		assertMatches(flipXY(expected), flipXY(getCliques(input)));
	}
	@Test
	public void single() {
		go(new TrapezoidGraphNode[] {N(1,2,3,4)}, N(1,2,3,4));
	}
	@Test
	public void no_overlap() {
		go(new TrapezoidGraphNode[] {N(1,2,1,2), N(1,2,3,4), N(3,4,1,2), N(3,4,3,4)}, N(1,2,1,2), N(1,2,3,4), N(3,4,1,2), N(3,4,3,4));
	}
	@Test
	public void contains() {
		go(new TrapezoidGraphNode[] {N(1,1,4,4), N(2,2,3,3)}, N(2,2,3,3,2));
	}
	@Test
	public void partial_overlap() {
		go(new TrapezoidGraphNode[] {N(1,1,4,4), N(0,2,2,3)}, N(1,2,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,1,4,4), N(0,2,2,5)}, N(1,2,2,4,2));
		go(new TrapezoidGraphNode[] {N(1,1,4,4), N(2,3,2,5)}, N(2,3,2,4,2));
		go(new TrapezoidGraphNode[] {N(1,1,4,4), N(3,5,2,5)}, N(3,4,2,4,2));
		go(new TrapezoidGraphNode[] {N(1,1,4,4), N(3,5,2,3)}, N(3,4,2,3,2));
	}
	@Test
	public void adjacent() {
		// adjacent X
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(5,6,0,5)}, N(1,4,1,4), N(5,6,0,5));
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(5,6,0,2)}, N(5,6,0,2), N(1,4,1,4));
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(5,6,2,3)}, N(5,6,2,3), N(1,4,1,4));
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(5,6,3,5)}, N(1,4,1,4), N(5,6,3,5));
		// adjacent Y
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(0,5,5,6)}, N(1,4,1,4), N(0,5,5,6));
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(0,2,5,6)}, N(1,4,1,4), N(0,2,5,6));
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(2,3,5,6)}, N(1,4,1,4), N(2,3,5,6));
		go(new TrapezoidGraphNode[] {N(1,4,1,4), N(3,5,5,6)}, N(1,4,1,4), N(3,5,5,6));
	}
	@Test
	public void x_start_overlap() {
		// only need to test first start1 <= second start1 due to symmetry 
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,1,1)}, N(1,1,1,1,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,1,2)}, N(1,1,1,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,1,3)}, N(1,1,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,1,4)}, N(1,1,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,2,2)}, N(1,1,2,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,2,3)}, N(1,1,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,2,4)}, N(1,1,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,3,3)}, N(1,1,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,3,4)}, N(1,1,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,1,4,4)});

		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,1,1)}, N(1,2,1,1,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,1,2)}, N(1,2,1,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,1,3)}, N(1,2,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,1,4)}, N(1,2,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,2,2)}, N(1,2,2,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,2,3)}, N(1,2,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,2,4)}, N(1,2,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,3,3)}, N(1,2,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,3,4)}, N(1,2,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,2,4,4)});
		
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,1,1)}, N(1,3,1,1,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,1,2)}, N(1,3,1,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,1,3)}, N(1,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,1,4)}, N(1,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,2,2)}, N(1,3,2,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,2,3)}, N(1,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,2,4)}, N(1,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,3,3)}, N(1,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,3,4)}, N(1,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,3,4,4)});
		
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,1,1)}, N(1,3,1,1,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,1,2)}, N(1,3,1,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,1,3)}, N(1,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,1,4)}, N(1,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,2,2)}, N(1,3,2,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,2,3)}, N(1,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,2,4)}, N(1,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,3,3)}, N(1,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,3,4)}, N(1,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(1,4,4,4)});
	}
	@Test
	public void x_end_overlap() {
		// first start1 < second start1
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,1,1)}, N(2,3,1,1,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,1,2)}, N(2,3,1,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,1,3)}, N(2,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,1,4)}, N(2,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,2,2)}, N(2,3,2,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,2,3)}, N(2,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,2,4)}, N(2,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,3,3)}, N(2,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,3,4)}, N(2,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(2,3,4,4)});
		
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,1,1)}, N(3,3,1,1,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,1,2)}, N(3,3,1,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,1,3)}, N(3,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,1,4)}, N(3,3,1,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,2,2)}, N(3,3,2,2,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,2,3)}, N(3,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,2,4)}, N(3,3,2,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,3,3)}, N(3,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,3,4)}, N(3,3,3,3,2));
		go(new TrapezoidGraphNode[] {N(1,3,1,3),N(3,3,4,4)});
	}
	@Test
	public void weighted() {
		go(new TrapezoidGraphNode[] {N(1,1,4,4,1), N(1,1,4,4,2)}, N(1,1,4,4,3));
	}
	@Test
	public void sequential_matching_weights() {
		go(new TrapezoidGraphNode[] { N(1, 2, 1, 1), N(3, 4, 1, 1), N(5, 6, 1, 1), N(7, 8, 1, 1)}, N(1, 2, 1, 1), N(3, 4, 1, 1), N(5, 6, 1, 1), N(7, 8, 1, 1));
	}
	@Test
	public void sequential_contained_cliques() {
		go(new TrapezoidGraphNode[] { N(1, 4, 1, 6), N(2, 3, 2, 3), N(2, 3, 4, 5) },
				N(2, 3, 2, 3, 2), N(2, 3, 4, 5, 2));
	}
	@Test
	public void clique_of_three() {
		go(new TrapezoidGraphNode[] { N(1,4,2,5), N(2,3,1,4), N(2,5,3,6) },
				N(2,3,3,4,3));
	}
	@Test
	public void many_cliques_order_n_squared() {
		// lattice pattern = maximal clique at every intersection
		//  | | | | | |  
		// -*-*-*-*-*-*-
		//  | | | | | |  
		// -*-*-*-*-*-*-
		//  | | | | | |  
		// -*-*-*-*-*-*-
		//  | | | | | |  
		// -*-*-*-*-*-*-
		//  | | | | | |  
		// -*-*-*-*-*-*-
		//  | | | | | |  
		// -*-*-*-*-*-*-
		//  | | | | | |  
		int n = 32; // technically n = # vertices = twice this number
		TrapezoidGraphNode[] expected = new TrapezoidGraphNode[n * n];
		TrapezoidGraphNode[] nodes = new TrapezoidGraphNode[2 * n];
		for (int i = 0; i < n; i++) {
			nodes[2 * i] = new TrapezoidGraphNode(2 * i, 2 * i, -1, 2 * n + 1, 1);
			nodes[2 * i + 1] = new TrapezoidGraphNode(-1, 2 * n + 1, 2 * i, 2 * i, 1);
			for (int j = 0; j < n; j++) {
				expected[n * i + j] = new TrapezoidGraphNode(2 * i, 2 * i, 2 * j, 2 * j, 2);
			}
		}
		go(nodes, expected);
	}
}
