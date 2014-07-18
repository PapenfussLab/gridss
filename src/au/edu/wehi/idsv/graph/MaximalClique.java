package au.edu.wehi.idsv.graph;

import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.SortedMap;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Queues;
/**
 * Rectangle graph implementation
 * 
 * D T Lee, Maximum clique problem of rectangle graphs, Advances in Computer Research, 1 (1983), pp. 91-107
 * 
 * @author Daniel Cameron
 */
public class MaximalClique {
	// Need a map since Set is missing a get method in Java and we don't want adding weight to an existing value to be an O(n) operation   
	private SortedMap<GraphNode, GraphNode> nodes = Maps.newTreeMap(GraphNode.ByStartXYEndXY); // sorted by startX, startY, then unimportant (but needed for uniqueness)
	/**
	 * Adds the given vertex to the trapezoid graph 
	 * @param node vertex to add
	 */
	public void add(GraphNode node) {
		GraphNode inMap = nodes.get(node);
		if (inMap == null) {
			inMap = new GraphNode(node.startX, node.endX, node.startY, node.endY, 0);
			nodes.put(inMap, inMap);
		}
		inMap.weight += node.weight;
	}
	private class MaximalCliqueEnumerator extends AbstractIterator<GraphNode> {
		private final PeekingIterator<GraphNode> nodeIt;
		private final Queue<GraphNode> outputBuffer = Queues.newArrayDeque();
		private final PriorityQueue<GraphNode> activeEndingX = new PriorityQueue<GraphNode>(11, GraphNode.ByEndXYStartXY); // sorted by endX
		private final CliqueScanlineActiveSet active = new CliqueScanlineActiveSet();
		public MaximalCliqueEnumerator(Iterator<GraphNode> nodeIt) {
			this.nodeIt = Iterators.peekingIterator(nodeIt);
		}
		/**
		 * Debugging checks asserting the internal data structures are correct 
		 */
		private void assertScanlineInvariants(long x) {
			if (nodeIt.hasNext() && nodeIt.peek().startX <= x) throw new RuntimeException(String.format("%s not yet processed by scanline %d", nodeIt.peek(), x));
			for (GraphNode n : nodes.keySet()) {
				boolean shouldBeActive = n.startX <= x && n.endX > x;
				if (shouldBeActive != activeEndingX.contains(n)) {
					throw new RuntimeException(String.format("%s unexpectedly %sin activeEndingX at %d (found %s)", n, shouldBeActive ? "not " : "", x, activeEndingX));
				}
				if (shouldBeActive) {
					active.assertContains(n);
				}
			}
		}
		public class CliqueScanlineActiveSet {
			private final NavigableMap<Long, YNode> active = Maps.newTreeMap();
			public CliqueScanlineActiveSet() {
				// sentinel node to reduce number of edge cases
				active.put(Long.MIN_VALUE, new YNode());
			}
			/**
			 * Debugging check asserting that the given node is in the active set as expected
			 * @param n
			 */
			public void assertContains(GraphNode n) {
				if (!active.containsKey(n.startY)) throw new RuntimeException(String.format("Missing y start node required by %s", n));
				if (!active.containsKey(n.endY + 1)) throw new RuntimeException(String.format("Missing y end node required by %s", n));
				if (active.get(n.startY).startHere == 0) throw new RuntimeException(String.format("Missing startHere>0 at y start node required by %s", n));
				if (active.floorEntry(n.endY).getValue().endHere == 0) throw new RuntimeException(String.format("Missing endHere>0 at y end node (%d) required by %s", active.floorEntry(n.endY).getKey() ,n));
			}
			private class YNode {
				public YNode() {
					startx = Long.MIN_VALUE;
					isMaximalClique = false;
					weight = 0;
					startHere = 0;
					endHere = 0;
				}
				public long startx;
				public boolean isMaximalClique;
				public float weight;
				public int startHere;
				public int endHere;
				@Override
				public String toString() {
					return String.format("startx=%d start#=%d end#=%d %s %f", startx, startHere, endHere, isMaximalClique ? "MAXIMAL" : "", weight);
				}
			}
			/**
			 * Adds a new node to the active set
			 * @param node
			 */
			public void add(GraphNode node) {
				split(node.startY);
				active.get(node.startY).startHere++;
				split(node.endY + 1);
				active.floorEntry(node.endY).getValue().endHere++;
				
				// Update the interval that our new nodes spans
				for (YNode ynode : active.subMap(node.startY, node.endY + 1).values()) {
					ynode.weight += node.weight;
					ynode.startx = node.startX;
					ynode.isMaximalClique = ynode.startHere > 0 && ynode.endHere > 0;
				}
			}
			/**
			 * Removes the given node from the active set, calling all maximal cliques this 
			 * node is part of that have not been called before
			 * @param node node to remove
			 * @return all new maximal cliques this node is part of
			 */
			public List<GraphNode> remove(GraphNode node) {
				// call maximal cliques
				List<GraphNode> cliques = Lists.newArrayList();
				for (PeekingIterator<Entry<Long, YNode>> it = Iterators.peekingIterator(active.subMap(node.startY, node.endY + 1).entrySet().iterator()); it.hasNext();){
					Entry<Long, YNode> entry = it.next();
					long starty = entry.getKey();
					YNode ynode = entry.getValue();
					if (ynode.isMaximalClique) {
						cliques.add(new GraphNode(ynode.startx, node.endX, starty, it.hasNext() ? it.peek().getKey() - 1 : node.endY, ynode.weight));
					}
					ynode.startx = Long.MIN_VALUE;
					ynode.isMaximalClique = false;
					// remove node contribution
					ynode.weight -= node.weight;
				}
				// Remove and coalese start/end bounds
				active.get(node.startY).startHere--;
				merge(node.startY);
				active.floorEntry(node.endY).getValue().endHere--;
				merge(node.endY + 1);
				return cliques;
			}
			/**
			 * Merges the given node with adjacent nodes (if possible)
			 * @param y
			 */
			private void merge(long y) {
				YNode node = active.get(y);
				if (node == null) throw new RuntimeException(String.format("Sanity check failure: node at position %d does not exist", y));
				YNode before = active.lowerEntry(y).getValue();
				if (before.endHere != 0 || node.startHere != 0) return; // Can't merge since they are still distinct
				if (before.isMaximalClique) throw new RuntimeException(String.format("Sanity check failure: partial maximal clique when merging %d", y));
				if (node.isMaximalClique) throw new RuntimeException(String.format("Sanity check failure: partial maximal clique when merging %d", y));
				if (before.weight != node.weight) throw new RuntimeException(String.format("Sanity check failure: nodes have different weights"));
				before.endHere = node.endHere;
				active.remove(y);
			}
			private void split(long y) {
				Entry<Long, YNode> entry = active.floorEntry(y);
				if (entry.getKey() == y) return; // Nothing to do - already split
				YNode before = entry.getValue();
				YNode node = new YNode();
				node.weight += before.weight;
				node.endHere = before.endHere;
				before.isMaximalClique = false;
				before.endHere = 0;
				active.put(y, node);
			}
			@Override
			public String toString() {
				StringBuilder sb = new StringBuilder();
				for (Entry<Long, YNode> entry : active.entrySet()) {
					sb.append(String.format("(starty=%d, %s)\n", entry.getValue(), entry.getKey()));
				}
				return sb.toString();
			}
		}
		/**
		 * Gets the next unprocessed scanline
		 * @return next scanline to process
		 */
		private long getNextX() {
			if (nodeIt.hasNext() && !activeEndingX.isEmpty()) {
				return Math.min(nodeIt.peek().startX, activeEndingX.peek().endX);
			} else if (nodeIt.hasNext()) {
				return nodeIt.peek().startX;
			} else if (!activeEndingX.isEmpty()) {
				return activeEndingX.peek().endX;
			} else {
				return Long.MAX_VALUE;
			}
		}
		private void processNextScanline() {
			long currentX = getNextX();
			// add new nodes to working set that start here
			while (nodeIt.hasNext() && nodeIt.peek().startX == currentX) {
				GraphNode startingHere = nodeIt.next();
				active.add(startingHere);
				activeEndingX.add(startingHere);
			}
			// remove nodes from working set that end here
			while (!activeEndingX.isEmpty() && activeEndingX.peek().endX == currentX) {
				GraphNode endingHere = activeEndingX.poll();
				outputBuffer.addAll(active.remove(endingHere));
			}
			assertScanlineInvariants(currentX);
		}
		@Override
		protected GraphNode computeNext() {
			while (outputBuffer.isEmpty() && getNextX() != Long.MAX_VALUE) {
				 processNextScanline();
			}
			if (!outputBuffer.isEmpty()) {
				return outputBuffer.poll();
			}
			return endOfData();
		}
	}
	/**
	 * Gets all maximal cliques in the graph
	 * @return maximal cliques ordered by endX, Y (startY and endY ordering are equivalent as maximal cliques cannot overlap)
	 */
	public Iterator<GraphNode> getAllMaximalCliques() {
		return new MaximalCliqueEnumerator(nodes.values().iterator());
	}
}
 