package au.edu.wehi.idsv.graph;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

public class MaximalCliqueTest {
	private GraphNode N(long startX, long endX, long startY, long endY) {
		return new GraphNode(startX, endX, startY, endY, 1);
	}
	private GraphNode N(long startX, long endX, long startY, long endY, int weight) {
		return new GraphNode(startX, endX, startY, endY, weight);
	}
	MaximalClique graph; 
	private GraphNode[] getCliques(GraphNode[] nodes) {
		graph = new MaximalClique();
		for (GraphNode n : nodes) {
			graph.add(n);
		}
		List<GraphNode> result = Lists.newArrayList(graph.getAllMaximalCliques());
		return result.toArray(new GraphNode[0]);
	}
	private GraphNode[] flipXY(GraphNode[] a) {
		GraphNode[] r = new GraphNode[a.length];
		for (int i = 0; i < a.length; i++) {
			r[i] = new GraphNode(a[i].startY, a[i].endY, a[i].startX, a[i].endX, a[i].evidence);
		}
		return r;
	}
	private String toString(GraphNode[] nodes) {
		StringBuilder sb = new StringBuilder();
		List<GraphNode> x = GraphNode.ByStartXYEndXY.sortedCopy(Arrays.asList(nodes));
		for (GraphNode n : x) {
			sb.append(n.toString());
			sb.append('\n');
		}
		return sb.toString();
	}
	private void assertMatches(GraphNode[] expected, GraphNode[] actual) {
		assertEquals(toString(expected), toString(actual));
	}
	private void go(GraphNode[] input, GraphNode... expected) {		
		assertMatches(expected, getCliques(input));
		assertMatches(flipXY(expected), getCliques(flipXY(input)));
	}
	@Test
	public void single() {
		go(new GraphNode[] {N(1,2,3,4,2)}, N(1,2,3,4,2));
	}
	@Test
	public void duplicate() {
		go(new GraphNode[] {N(1,2,3,4), N(1,2,3,4)}, N(1,2,3,4,2));
	}
	@Test
	public void no_overlap() {
		go(new GraphNode[] {N(1,2,1,2), N(1,2,3,4), N(3,4,1,2), N(3,4,3,4)}, N(1,2,1,2), N(1,2,3,4), N(3,4,1,2), N(3,4,3,4));
	}
	@Test
	public void contains() {
		go(new GraphNode[] {N(1,4,1,4), N(2,2,3,3)}, N(2,2,3,3,2));
	}
	@Test
	public void contains_subsequent() {
		// smaller
		go(new GraphNode[] {N(1,6,1,6), N(2,5,2,3), N(3,4,4,5)}, N(2,5,2,3,2), N(3,4,4,5,2));
		// larger
		go(new GraphNode[] {N(1,6,1,6), N(2,5,4,5), N(3,4,2,3)}, N(2,5,4,5,2), N(3,4,2,3,2));
		// overlap
		go(new GraphNode[] {N(1,6,1,6), N(2,4,4,5), N(3,5,2,3)}, N(2,4,4,5,2), N(3,5,2,3,2));
		go(new GraphNode[] {N(1,6,1,6), N(3,5,4,5), N(2,4,2,3)}, N(3,5,4,5,2), N(2,4,2,3,2));
	}
	@Test
	public void adjacent_parents() {
		go(new GraphNode[] {N(1,3,1,5), N(3,5,1,5), N(2,4,1,2), N(2,4,3,4)}, N(3,3,1,2,3), N(3,3,3,4,3));
		go(new GraphNode[] {N(1,3,1,6), N(4,6,1,6), N(2,5,2,3), N(2,5,4,5)}, N(2,3,2,3,2), N(4,5,2,3,2), N(2,3,4,5,2), N(4,5,4,5,2));
	}
	@Test
	public void partial_overlap() {
		go(new GraphNode[] {N(1,4,1,4), N(0,2,2,3)}, N(1,2,2,3,2));
		go(new GraphNode[] {N(1,4,1,4), N(0,2,2,5)}, N(1,2,2,4,2));
		go(new GraphNode[] {N(1,4,1,4), N(2,3,2,5)}, N(2,3,2,4,2));
		go(new GraphNode[] {N(1,4,1,4), N(3,5,2,5)}, N(3,4,2,4,2));
		go(new GraphNode[] {N(1,4,1,4), N(3,5,2,3)}, N(3,4,2,3,2));
	}
	@Test
	public void adjacent() {
		// adjacent X
		go(new GraphNode[] {N(1,4,1,4), N(5,6,0,5)}, N(1,4,1,4), N(5,6,0,5));
		go(new GraphNode[] {N(1,4,1,4), N(5,6,0,2)}, N(5,6,0,2), N(1,4,1,4));
		go(new GraphNode[] {N(1,4,1,4), N(5,6,2,3)}, N(5,6,2,3), N(1,4,1,4));
		go(new GraphNode[] {N(1,4,1,4), N(5,6,3,5)}, N(1,4,1,4), N(5,6,3,5));
		// adjacent Y
		go(new GraphNode[] {N(1,4,1,4), N(0,5,5,6)}, N(1,4,1,4), N(0,5,5,6));
		go(new GraphNode[] {N(1,4,1,4), N(0,2,5,6)}, N(1,4,1,4), N(0,2,5,6));
		go(new GraphNode[] {N(1,4,1,4), N(2,3,5,6)}, N(1,4,1,4), N(2,3,5,6));
		go(new GraphNode[] {N(1,4,1,4), N(3,5,5,6)}, N(1,4,1,4), N(3,5,5,6));
	}
	@Test
	public void start_overlap() {
		// only need to test first startX <= second startX due to symmetry 
		go(new GraphNode[] {N(1,3,1,3),N(1,1,1,1)}, N(1,1,1,1,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,1,2)}, N(1,1,1,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,1,3)}, N(1,1,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,1,4)}, N(1,1,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,2,2)}, N(1,1,2,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,2,3)}, N(1,1,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,2,4)}, N(1,1,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,3,3)}, N(1,1,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,3,4)}, N(1,1,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,1,4,4)}, N(1,3,1,3),N(1,1,4,4));

		go(new GraphNode[] {N(1,3,1,3),N(1,2,1,1)}, N(1,2,1,1,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,1,2)}, N(1,2,1,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,1,3)}, N(1,2,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,1,4)}, N(1,2,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,2,2)}, N(1,2,2,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,2,3)}, N(1,2,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,2,4)}, N(1,2,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,3,3)}, N(1,2,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,3,4)}, N(1,2,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,2,4,4)}, N(1,3,1,3),N(1,2,4,4));
		
		go(new GraphNode[] {N(1,3,1,3),N(1,3,1,1)}, N(1,3,1,1,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,1,2)}, N(1,3,1,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,1,3)}, N(1,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,1,4)}, N(1,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,2,2)}, N(1,3,2,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,2,3)}, N(1,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,2,4)}, N(1,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,3,3)}, N(1,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,3,4)}, N(1,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,3,4,4)}, N(1,3,1,3),N(1,3,4,4));
		
		go(new GraphNode[] {N(1,3,1,3),N(1,4,1,1)}, N(1,3,1,1,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,1,2)}, N(1,3,1,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,1,3)}, N(1,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,1,4)}, N(1,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,2,2)}, N(1,3,2,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,2,3)}, N(1,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,2,4)}, N(1,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,3,3)}, N(1,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,3,4)}, N(1,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(1,4,4,4)}, N(1,3,1,3),N(1,4,4,4));
	}
	@Test
	public void end_overlap() {
		// first startX < second startX
		go(new GraphNode[] {N(1,3,1,3),N(2,3,1,1)}, N(2,3,1,1,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,1,2)}, N(2,3,1,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,1,3)}, N(2,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,1,4)}, N(2,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,2,2)}, N(2,3,2,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,2,3)}, N(2,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,2,4)}, N(2,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,3,3)}, N(2,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,3,4)}, N(2,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(2,3,4,4)}, N(1,3,1,3),N(2,3,4,4));
		
		go(new GraphNode[] {N(1,3,1,3),N(3,3,1,1)}, N(3,3,1,1,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,1,2)}, N(3,3,1,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,1,3)}, N(3,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,1,4)}, N(3,3,1,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,2,2)}, N(3,3,2,2,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,2,3)}, N(3,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,2,4)}, N(3,3,2,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,3,3)}, N(3,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,3,4)}, N(3,3,3,3,2));
		go(new GraphNode[] {N(1,3,1,3),N(3,3,4,4)}, N(1,3,1,3),N(3,3,4,4));
	}
	@Test
	public void weighted() {
		go(new GraphNode[] {N(4,4,1,1,1), N(4,4,1,1,2)}, N(4,4,1,1,3));
		go(new GraphNode[] {N(1,1,4,4,1), N(1,1,4,4,2)}, N(1,1,4,4,3));
	}
	@Test
	public void sequential_matching_weights() {
		go(new GraphNode[] { N(1, 2, 1, 1), N(3, 4, 1, 1), N(5, 6, 1, 1), N(7, 8, 1, 1)}, N(1, 2, 1, 1), N(3, 4, 1, 1), N(5, 6, 1, 1), N(7, 8, 1, 1));
	}
	@Test
	public void sequential_contained_cliques() {
		go(new GraphNode[] { N(1, 4, 1, 6), N(2, 3, 2, 3), N(2, 3, 4, 5) },
				N(2, 3, 2, 3, 2), N(2, 3, 4, 5, 2));
	}
	@Test
	public void clique_of_three() {
		go(new GraphNode[] { N(1,4,2,5), N(2,3,1,4), N(2,5,3,6) },
				N(2,3,3,4,3));
	}
	/**
	 * lattice pattern = maximal clique at every intersection
	 *  | | | | | |  
	 * -*-*-*-*-*-*-
	 *  | | | | | |  
	 * -*-*-*-*-*-*-
	 *  | | | | | |  
	 * -*-*-*-*-*-*-
	 *  | | | | | |  
	 * -*-*-*-*-*-*-
	 *  | | | | | |  
	 * -*-*-*-*-*-*-
	 *  | | | | | |  
	 * -*-*-*-*-*-*-
	 *  | | | | | |  
	 * @param n number of lattice rows and columns
	 */
	private void test_lattice(int n) {
		GraphNode[] expected = new GraphNode[n * n];
		GraphNode[] nodes = new GraphNode[2 * n];
		for (int i = 0; i < n; i++) {
			nodes[2 * i] = N(2 * i, 2 * i, -1, 2 * n + 1, 1);
			nodes[2 * i + 1] = N(-1, 2 * n + 1, 2 * i, 2 * i, 1);
			for (int j = 0; j < n; j++) {
				expected[n * i + j] = N(2 * i, 2 * i, 2 * j, 2 * j, 2);
			}
		}
		go(nodes, expected);
	}
	@Test
	public void many_cliques_order_n_squared() {
		for (int i = 0; i < 8; i++) {
			test_lattice(1 << i);			
		}
	}
}
