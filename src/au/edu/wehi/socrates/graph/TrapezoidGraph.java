package au.edu.wehi.socrates.graph;

import java.util.Iterator;
import java.util.List;
import java.util.NavigableMap;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.SortedSet;

import javafx.collections.transformation.SortedList;
import au.edu.wehi.socrates.BreakpointInterval;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Queues;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
/**
 * Trapezoid graph implementation
 * 
 * @author Daniel Cameron
 */
public class TrapezoidGraph {
	//private SortedSetMultimap<Long, TrapezoidGraphNode> intervals = TreeMultimap.create(null, null);
	private SortedSet<TrapezoidGraphNode> nodes; // sorted by (start1, start2)
	/**
	 * Adds the given vertex to the trapezoid graph 
	 * @param node vertex to add
	 */
	public void add(TrapezoidGraphNode node) {
		if (nodes.contains(node)) {
			// combine weights
			throw new RuntimeException("NYI");
		} else {
			nodes.add(node);
		}
	}
	private class MaximalCliqueEnumerator extends AbstractIterator<TrapezoidGraphNode> {
		private PeekingIterator<TrapezoidGraphNode> nodeIt;
		private final Queue<TrapezoidGraphNode> outputBuffer = Queues.newArrayDeque();
		PriorityQueue<TrapezoidGraphNode> nodeEnding; // sorted by (end1, start2)
		//private SortedList<SortedList<TrapezoidGraphNode>> nodeEnding; // sorted by (end1); sorted by (start2)
		private SortedSet<TrapezoidGraphNode> activeNodes; // sorted by (start2)
		/**
		 * Gets the next unprocessed scanline
		 * @return next scanline to process
		 */
		private long getNextScanline() {
			if (nodeIt.hasNext() && !nodeEnding.isEmpty()) {
				return Math.min(nodeIt.peek().start1, nodeEnding.peek().end1);
			} else if (nodeIt.hasNext()) {
				return nodeIt.peek().start1;
			} else if (!nodeEnding.isEmpty()) {
				return nodeEnding.peek().end1;
			} else {
				return Long.MAX_VALUE;
			}
		}
		private void addToActive(TrapezoidGraphNode node) {
			activeNodes.add(node);
			nodeEnding.add(node);
		}
		private void processNextScanline() {
			long currentLine = getNextScanline();
			// add new nodes to working set
				// nodeEndings += new
				// activeNodes += new
			while (nodeIt.hasNext() && nodeIt.peek().start1 == currentLine) {
				addToActive(nodeIt.next());
			}			
			List<TrapezoidGraphNode> endingHere = pollScanlineEnd1(nodeEnding, currentLine);
			// IntervalSet<Long, Long> callIntervals = getScanlines(endingHere);
			// for maximal cliques in activeNodes intersect callIntervals
				//outputBuffer.add(call)
			
			// for node in endingHere
				// remove node from activeNodes
		}
		@Override
		protected TrapezoidGraphNode computeNext() {
			while (!outputBuffer.isEmpty()) {
				 processNextScanline();
			 }
			if (outputBuffer.isEmpty()) {
				return endOfData();
			} else {
				return outputBuffer.poll();
			}
		}
	}
	/**
	 * removes items from the head of the queue whose end1 matches the given value
	 * @param queue queue to poll
	 * @param end1 end1 to match 
	 * @return list of all matches removed from the queue
	 */
	private static List<TrapezoidGraphNode> pollScanlineEnd1(Queue<TrapezoidGraphNode> queue, long end1) {
		List<TrapezoidGraphNode> list = Lists.newArrayList();
		while (queue.peek().end1 == end1) {
			list.add(queue.poll());
		}
		return list;
	}
	public Iterator<TrapezoidGraphNode> getAllMaximalCliques() {
		return new MaximalCliqueEnumerator();
		// Traverse graph and find all maximum cliques
		// Boxicity = 2
		
		// geometric equivalence:
				// maximal clique == area (ie all points (p1, p2)) in trapezoid of overlap 
				// maximal clique equivalent to all reads matching: a given point
				
				// OEA evidence interval for fusion requires mate to overlap actual fusion sequence [minFragmentSize, maxFragmentSize]
				// OEA evidence interval is [0, maxFragmentSize]
	}
}
 