package au.edu.wehi.idsv.graph;

import au.edu.wehi.idsv.graph.MaximumCliqueIntervalGraph.Node;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;


public class MaximumCliqueIntervalGraphTest {
	@Test
	public void should_return_lowest_position_best_score() {
		List<Node> nodes = new ArrayList<Node>();
		nodes.add(new Node(1, 2, 1));
		nodes.add(new Node(2, 3, 1));
		nodes.add(new Node(3, 4, 1));
		nodes.add(new Node(4, 5, 1));
		nodes.add(new Node(5, 6, 1));
		Node result = new MaximumCliqueIntervalGraph().calculateMaximumClique(nodes);
		assertEquals(2L, result.start);
		assertEquals(2L, result.stop);
		assertEquals(2L, result.weight);
	}
	@Test
	public void should_incorporate_node_weight() {
		List<Node> nodes = new ArrayList<Node>();
		nodes.add(new Node(1, 2, 3));
		nodes.add(new Node(3, 4, 1));
		nodes.add(new Node(3, 4, 1));
		Node result = new MaximumCliqueIntervalGraph().calculateMaximumClique(nodes); 
		assertEquals(1L, result.start);
		assertEquals(2L, result.stop);
		assertEquals(3L, result.weight);
	}
	@Test
	public void should_return_best_scoring_clique() {
		List<Node> nodes = new ArrayList<Node>();
		nodes.add(new Node(1, 2, 1));
		nodes.add(new Node(2, 3, 1));
		nodes.add(new Node(3, 4, 2));
		nodes.add(new Node(4, 5, 1));
		nodes.add(new Node(5, 6, 1));
		nodes.add(new Node(4, 10, 1));
		Node result = new MaximumCliqueIntervalGraph().calculateMaximumClique(nodes);
		assertEquals(4L, result.start);
		assertEquals(4L, result.stop);
		assertEquals(4L, result.weight);
		
		nodes = new ArrayList<Node>();
		nodes.add(new Node(10, 15, 1));
		nodes.add(new Node(15, 18, 1));
		nodes.add(new Node(17, 18, 1));
		nodes.add(new Node(9, 15, 1));
		nodes.add(new Node(19, 20, 1));
		nodes.add(new Node(19, 21, 1));
		result = new MaximumCliqueIntervalGraph().calculateMaximumClique(nodes);
		assertEquals(15, result.start);
		assertEquals(15, result.stop);
		assertEquals(3, result.weight);
	}
}
