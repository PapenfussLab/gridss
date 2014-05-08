package au.edu.wehi.socrates.graph;

import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.SortedMap;
import java.util.SortedSet;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Queues;
import com.google.common.collect.Sets;
/**
 * Trapezoid graph implementation
 * 
 * @author Daniel Cameron
 */
public class TrapezoidGraph {
	// Need a map since Set is missing a get method in Java and we don't want adding weight to an existing value to be an O(n) operation   
	private SortedMap<TrapezoidGraphNode, TrapezoidGraphNode> nodes = Maps.newTreeMap(new TrapezoidGraphNodeStartXYComparator()); // sorted by startX, startY, then unimportant (but needed for uniqueness)
	/**
	 * Adds the given vertex to the trapezoid graph 
	 * @param node vertex to add
	 */
	public void add(TrapezoidGraphNode node) {
		if (nodes.containsKey(node)) {
			node = new TrapezoidGraphNode(node, nodes.get(node).weight);
		}
		nodes.put(node, node);
	}
	private static class MaximalCliqueEnumerator extends AbstractIterator<TrapezoidGraphNode> {
		private final PeekingIterator<TrapezoidGraphNode> nodeIt;
		private final Queue<TrapezoidGraphNode> outputBuffer = Queues.newArrayDeque();
		private final PriorityQueue<TrapezoidGraphNode> activeEndingX = new PriorityQueue<TrapezoidGraphNode>(11, new TrapezoidGraphNodeEndXYComparator()); // sorted by endX
		private final SortedSet<TrapezoidGraphNode> activeX = Sets.newTreeSet(new TrapezoidGraphNodeStartYXComparator()); // sorted by startY, then unimportant (but needed for uniqueness)
		public MaximalCliqueEnumerator(Iterator<TrapezoidGraphNode> nodeIt) {
			this.nodeIt = Iterators.peekingIterator(nodeIt);
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
			// add new nodes to working set
			while (nodeIt.hasNext() && nodeIt.peek().startX == currentX) {
				TrapezoidGraphNode startingHere = activeEndingX.poll();
				activeX.add(startingHere);
				activeEndingX.add(startingHere);
			}
			// process maximal cliques on current scan-line
			outputBuffer.addAll(maximalCliques(currentX));
			// remove cliques
			while (!activeEndingX.isEmpty() && activeEndingX.peek().endX == currentX) {
				TrapezoidGraphNode endingHere = activeEndingX.poll();
				activeX.remove(endingHere);
			}
		}
		@Override
		protected TrapezoidGraphNode computeNext() {
			while (outputBuffer.isEmpty() && getNextX() != Long.MAX_VALUE) {
				 processNextScanline();
			}
			if (outputBuffer.isEmpty()) {
				return endOfData();
			}
			return outputBuffer.poll();
		}
		private List<TrapezoidGraphNode> maximalCliques(long xend) {
			List<TrapezoidGraphNode> cliques = Lists.newArrayList();
			SortedSet<TrapezoidGraphNode> xstart = Sets.newTreeSet(new TrapezoidGraphNodeStartXYComparator()); // sorted by startX, then unimportant (but needed for uniqueness)
			PriorityQueue<TrapezoidGraphNode> endY = new PriorityQueue<TrapezoidGraphNode>(11, new TrapezoidGraphNodeEndYXComparator()); // sorted by endY
			long ystart = Long.MIN_VALUE;
			boolean inMaximal = false;
			double currentWeight = 0;
			// TODO: can we move the calling into the endingX loop of processNextScanline?
			// what data structure would be need? Some sort of hierarchical range interval set
			PeekingIterator<TrapezoidGraphNode> yit = Iterators.peekingIterator(activeX.iterator());
			for (long y = getNextY(yit, endY); y != Long.MAX_VALUE; y = getNextY(yit, endY)) {
				// add starting set
				while (yit.hasNext() && yit.peek().startY == y) {
					TrapezoidGraphNode startingHere = yit.next();
					endY.add(startingHere);
					xstart.add(startingHere);
					currentWeight += startingHere.weight;
					if (startingHere.endX == xend) {
						// upper x bound of maximal clique must match upper x bound of a clique member
						// only set maximal if 
						inMaximal = true;
						ystart = y;
					}
				}
				if (endY.peek().endY == y) {
					if (inMaximal) {
						// we just hit the end of a maximal clique
						cliques.add(new TrapezoidGraphNode(xstart.first().startX, xend, ystart, y, currentWeight));
					}
					// remove all nodes that end here from the active set
					while (endY.peek().endY == y) {
						TrapezoidGraphNode endingHere = endY.poll();
						xstart.remove(endingHere);
						currentWeight -= endingHere.weight;
					}
					inMaximal = false; // no longer maximal
					ystart = Long.MIN_VALUE;
				}
			}
			return cliques;
		}
	}
	private static long getNextY(PeekingIterator<TrapezoidGraphNode> startIterator, Queue<TrapezoidGraphNode> endCollection) {
		if (startIterator.hasNext() && !endCollection.isEmpty()) {
			return Math.min(startIterator.peek().startY, endCollection.peek().endY);
		} else if (startIterator.hasNext()) {
			return startIterator.peek().startY;
		} else if (!endCollection.isEmpty()) {
			return endCollection.peek().endY;
		} else {
			return Long.MAX_VALUE;
		}
	}
	/**
	 * Gets all maximal cliques in the graph
	 * @return maximal cliques
	 */
	public Iterator<TrapezoidGraphNode> getAllMaximalCliques() {
		return new MaximalCliqueEnumerator(nodes.values().iterator());
	}
}
 