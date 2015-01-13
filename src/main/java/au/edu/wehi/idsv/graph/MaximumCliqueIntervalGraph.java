package au.edu.wehi.idsv.graph;

import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

/**
 * Calculates the maximum clique of an interval graph
 * @author cameron.d
 *
 */
public class MaximumCliqueIntervalGraph {
	public static class Node {
		public final long start;
		public final long stop;
		public final long weight;
		/**
		 * Interval graph node
		 * @param start start position, inclusive
		 * @param stop end position, inclusive
		 * @param weight node weight
		 */
		public Node(long start, long stop, long weight) {
			assert(weight > 0);
			this.start = start;
			this.stop = stop;
			this.weight = weight;
		}
		private Node() {
			this.start = 0;
			this.stop = 0;
			this.weight = 0;
		}
		@Override
		public String toString() {
			return String.format("(%d, %d, %d)", start, stop, weight);
		}
	}
	private static final Ordering<Node> ByStop = new Ordering<Node>() {
		public int compare(Node arg1, Node arg2) {
			return Longs.compare(arg1.stop, arg2.stop);
		}
	};
	private static final Ordering<Node> ByStart = new Ordering<Node>() {
		public int compare(Node arg1, Node arg2) {
			return Longs.compare(arg1.start, arg2.start);
		}
	};
	private PriorityQueue<Node> active = new PriorityQueue<Node>(16, ByStop);
	private Node best = new Node();
	private long activeStart = Long.MIN_VALUE;
	private long activeWeight = 0;
	/**
	 * Calculates the maximum clique for the given interval graph
	 * @param nodes weighted intervals
	 * @return maximum clique
	 */
	public Node calculateMaximumClique(List<Node> nodes) {
		Collections.sort(nodes, ByStart);
		for (Node startNode : nodes) {
			processStops(startNode.start);
			// process start node
			activeStart = startNode.start;
			activeWeight += startNode.weight;
			active.add(startNode);
		}
		processStops(Long.MAX_VALUE);
		return best;
	}
	private void processStops(long before) {
		while (!active.isEmpty() && active.peek().stop < before) {
			Node stopNode = active.poll();
			long stop = stopNode.stop;
			long stopWeight = stopNode.weight;
			while (!active.isEmpty() && active.peek().stop == stop) {
				stopWeight += active.poll().weight;
			}
			Node clique = new Node(activeStart, stop, activeWeight);
			addMaximalClique(clique);
			activeStart = stop + 1;
			activeWeight -= stopWeight;
		}
	}
	private void addMaximalClique(Node clique) {
		if (clique.weight > best.weight) {
			best = clique;
		}
	}
}
