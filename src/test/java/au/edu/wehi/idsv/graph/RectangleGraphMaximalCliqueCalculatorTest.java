package au.edu.wehi.idsv.graph;

import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class RectangleGraphMaximalCliqueCalculatorTest {
	private RectangleGraphNode N(long startX, long endX, long startY, long endY) {
		return new RectangleGraphNode(startX, endX, startY, endY, 1, 1);
	}
	private RectangleGraphNode N(long startX, long endX, long startY, long endY, int weight, int exactWeight) {
		return new RectangleGraphNode(startX, endX, startY, endY, weight, exactWeight);
	}
	private RectangleGraphNode N(long startX, long endX, long startY, long endY, int weight) {
		return new RectangleGraphNode(startX, endX, startY, endY, weight, weight);
	}
	RectangleGraphMaximalCliqueCalculator graph; 
	private RectangleGraphNode[] getCliques(RectangleGraphNode[] nodes) {
		Arrays.sort(nodes, 0, nodes.length, RectangleGraphNode.ByStartXYEndXY);
		graph = new RectangleGraphMaximalCliqueCalculator();
		List<RectangleGraphNode> result = Lists.newArrayList();
		for (int i = 0; i < nodes.length; i++) {
			result.addAll(graph.next(nodes[i]));
		}
		result.addAll(graph.complete());
		return result.toArray(new RectangleGraphNode[0]);
	}
	private RectangleGraphNode[] flipXY(RectangleGraphNode[] a) {
		RectangleGraphNode[] r = new RectangleGraphNode[a.length];
		for (int i = 0; i < a.length; i++) {
			r[i] = new RectangleGraphNode(a[i].startY, a[i].endY, a[i].startX, a[i].endX, a[i].weight, a[i].exactWeight);
		}
		return r;
	}
	private String toString(RectangleGraphNode[] nodes) {
		StringBuilder sb = new StringBuilder();
		List<RectangleGraphNode> x = RectangleGraphNode.ByStartXYEndXY.sortedCopy(Arrays.asList(nodes));
		for (RectangleGraphNode n : x) {
			sb.append(n.toString());
			sb.append('\n');
		}
		return sb.toString();
	}
	private void assertMatches(RectangleGraphNode[] expected, RectangleGraphNode[] actual) {
		assertEquals(toString(expected), toString(actual));
	}
	private void go(RectangleGraphNode[] input, RectangleGraphNode... expected) {		
		assertMatches(expected, getCliques(input));
		assertMatches(flipXY(expected), getCliques(flipXY(input)));
	}
	@Test
	public void single_base() {
		go(new RectangleGraphNode[] {N(1,1,1,1,2)}, N(1,1,1,1,2));
	}
	@Test
	public void single_base_adj() {
		go(new RectangleGraphNode[] {N(1,1,1,1), N(1,1,2,2)}, N(1,1,1,1), N(1,1,2,2));
	}
	@Test
	public void single_base_overlap() {
		go(new RectangleGraphNode[] {N(1,1,1,1), N(1,1,1,2)}, N(1,1,1,1,2));
	}
	@Test
	public void single() {
		go(new RectangleGraphNode[] {N(1,2,3,4,2)}, N(1,2,3,4,2));
	}
	@Test
	public void duplicate() {
		go(new RectangleGraphNode[] {N(1,2,3,4), N(1,2,3,4)}, N(1,2,3,4,2));
	}
	@Test
	public void no_overlap() {
		go(new RectangleGraphNode[] {N(1,2,1,2), N(1,2,3,4), N(3,4,1,2), N(3,4,3,4)}, N(1,2,1,2), N(1,2,3,4), N(3,4,1,2), N(3,4,3,4));
	}
	@Test
	public void contains() {
		go(new RectangleGraphNode[] {N(1,4,1,4), N(2,2,3,3)}, N(2,2,3,3,2));
	}
	@Test
	public void contains_subsequent() {
		// smaller
		go(new RectangleGraphNode[] {N(1,6,1,6), N(2,5,2,3), N(3,4,4,5)}, N(2,5,2,3,2), N(3,4,4,5,2));
		// larger
		go(new RectangleGraphNode[] {N(1,6,1,6), N(2,5,4,5), N(3,4,2,3)}, N(2,5,4,5,2), N(3,4,2,3,2));
		// overlap
		go(new RectangleGraphNode[] {N(1,6,1,6), N(2,4,4,5), N(3,5,2,3)}, N(2,4,4,5,2), N(3,5,2,3,2));
		go(new RectangleGraphNode[] {N(1,6,1,6), N(3,5,4,5), N(2,4,2,3)}, N(3,5,4,5,2), N(2,4,2,3,2));
	}
	@Test
	public void adjacent_parents() {
		go(new RectangleGraphNode[] {N(1,3,1,5), N(3,5,1,5), N(2,4,1,2), N(2,4,3,4)}, N(3,3,1,2,3), N(3,3,3,4,3));
		go(new RectangleGraphNode[] {N(1,3,1,6), N(4,6,1,6), N(2,5,2,3), N(2,5,4,5)}, N(2,3,2,3,2), N(4,5,2,3,2), N(2,3,4,5,2), N(4,5,4,5,2));
	}
	@Test
	public void partial_overlap() {
		go(new RectangleGraphNode[] {N(1,4,1,4), N(0,2,2,3)}, N(1,2,2,3,2));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(0,2,2,5)}, N(1,2,2,4,2));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(2,3,2,5)}, N(2,3,2,4,2));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(3,5,2,5)}, N(3,4,2,4,2));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(3,5,2,3)}, N(3,4,2,3,2));
	}
	@Test
	public void adjacent() {
		// adjacent X
		go(new RectangleGraphNode[] {N(1,4,1,4), N(5,6,0,5)}, N(1,4,1,4), N(5,6,0,5));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(5,6,0,2)}, N(5,6,0,2), N(1,4,1,4));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(5,6,2,3)}, N(5,6,2,3), N(1,4,1,4));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(5,6,3,5)}, N(1,4,1,4), N(5,6,3,5));
		// adjacent Y
		go(new RectangleGraphNode[] {N(1,4,1,4), N(0,5,5,6)}, N(1,4,1,4), N(0,5,5,6));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(0,2,5,6)}, N(1,4,1,4), N(0,2,5,6));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(2,3,5,6)}, N(1,4,1,4), N(2,3,5,6));
		go(new RectangleGraphNode[] {N(1,4,1,4), N(3,5,5,6)}, N(1,4,1,4), N(3,5,5,6));
	}
	@Test
	public void start_overlap() {
		// only need to test first startX <= second startX due to symmetry 
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,1,1)}, N(1,1,1,1,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,1,2)}, N(1,1,1,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,1,3)}, N(1,1,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,1,4)}, N(1,1,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,2,2)}, N(1,1,2,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,2,3)}, N(1,1,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,2,4)}, N(1,1,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,3,3)}, N(1,1,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,3,4)}, N(1,1,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,1,4,4)}, N(1,3,1,3),N(1,1,4,4));

		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,1,1)}, N(1,2,1,1,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,1,2)}, N(1,2,1,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,1,3)}, N(1,2,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,1,4)}, N(1,2,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,2,2)}, N(1,2,2,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,2,3)}, N(1,2,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,2,4)}, N(1,2,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,3,3)}, N(1,2,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,3,4)}, N(1,2,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,2,4,4)}, N(1,3,1,3),N(1,2,4,4));
		
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,1,1)}, N(1,3,1,1,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,1,2)}, N(1,3,1,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,1,3)}, N(1,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,1,4)}, N(1,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,2,2)}, N(1,3,2,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,2,3)}, N(1,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,2,4)}, N(1,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,3,3)}, N(1,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,3,4)}, N(1,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,3,4,4)}, N(1,3,1,3),N(1,3,4,4));
		
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,1,1)}, N(1,3,1,1,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,1,2)}, N(1,3,1,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,1,3)}, N(1,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,1,4)}, N(1,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,2,2)}, N(1,3,2,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,2,3)}, N(1,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,2,4)}, N(1,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,3,3)}, N(1,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,3,4)}, N(1,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(1,4,4,4)}, N(1,3,1,3),N(1,4,4,4));
	}
	@Test
	public void end_overlap() {
		// first startX < second startX
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,1,1)}, N(2,3,1,1,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,1,2)}, N(2,3,1,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,1,3)}, N(2,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,1,4)}, N(2,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,2,2)}, N(2,3,2,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,2,3)}, N(2,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,2,4)}, N(2,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,3,3)}, N(2,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,3,4)}, N(2,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(2,3,4,4)}, N(1,3,1,3),N(2,3,4,4));
		
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,1,1)}, N(3,3,1,1,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,1,2)}, N(3,3,1,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,1,3)}, N(3,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,1,4)}, N(3,3,1,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,2,2)}, N(3,3,2,2,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,2,3)}, N(3,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,2,4)}, N(3,3,2,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,3,3)}, N(3,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,3,4)}, N(3,3,3,3,2));
		go(new RectangleGraphNode[] {N(1,3,1,3),N(3,3,4,4)}, N(1,3,1,3),N(3,3,4,4));
	}
	@Test
	public void weighted() {
		go(new RectangleGraphNode[] {N(4,4,1,1,1), N(4,4,1,1,2)}, N(4,4,1,1,3));
		go(new RectangleGraphNode[] {N(1,1,4,4,1), N(1,1,4,4,2)}, N(1,1,4,4,3));
	}
	@Test
	public void sequential_matching_weights() {
		go(new RectangleGraphNode[] { N(1, 2, 1, 1), N(3, 4, 1, 1), N(5, 6, 1, 1), N(7, 8, 1, 1)}, N(1, 2, 1, 1), N(3, 4, 1, 1), N(5, 6, 1, 1), N(7, 8, 1, 1));
	}
	@Test
	public void sequential_contained_cliques() {
		go(new RectangleGraphNode[] { N(1, 4, 1, 6), N(2, 3, 2, 3), N(2, 3, 4, 5) },
				N(2, 3, 2, 3, 2), N(2, 3, 4, 5, 2));
	}
	@Test
	public void clique_of_three() {
		go(new RectangleGraphNode[] { N(1,4,2,5), N(2,3,1,4), N(2,5,3,6) },
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
		RectangleGraphNode[] expected = new RectangleGraphNode[n * n];
		RectangleGraphNode[] nodes = new RectangleGraphNode[2 * n];
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
	@Test
	public void exhaustive() {
		int size = 8;
		List<RectangleGraphNode> nodes = new ArrayList<RectangleGraphNode>();
		for (int startx = 0; startx < size; startx++) {
			for (int endx = startx; endx < size; endx++) {
				for (int starty = 0; starty < size; starty++) {
					for (int endy = starty; endy < size; endy++) {
						nodes.add(new RectangleGraphNode(startx, endx, starty, endy, 1, 0));
					}
				}
			}
		}
		RectangleGraphNode[] cliques = getCliques(nodes.toArray(new RectangleGraphNode[nodes.size()]));
		assertEquals(size * size, cliques.length); // clique at every grid position
	}
	@Test
	public void should_handle_greater_than_int_max_value() {
		int size = 3;
		List<RectangleGraphNode> nodes = new ArrayList<RectangleGraphNode>();
		for (int startx = 0; startx < size; startx++) {
			for (int endx = startx; endx < size; endx++) {
				for (int starty = 0; starty < size; starty++) {
					for (int endy = starty; endy < size; endy++) {
						nodes.add(new RectangleGraphNode(startx, endx, starty, endy, Integer.MAX_VALUE - 10, 0));
					}
				}
			}
		}
		RectangleGraphNode[] cliques = getCliques(nodes.toArray(new RectangleGraphNode[nodes.size()]));
		assertEquals(size * size, cliques.length); // clique at every grid position
	}
	@Test
	public void should_treat_weight_and_exactWeight_as_two_separate_weightings() {
		go(new RectangleGraphNode[]{
						N(1, 4, 1, 1, 10, 5),
						N(2, 3, 1, 1, 5, 2),
						N(4, 4, 1, 1, 1, 1)
				},
				N(2, 3, 1, 1, 15, 7),
				N(4, 4, 1, 1, 11, 6));
	}
}
