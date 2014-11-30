package au.edu.wehi.idsv.graph;

import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.Defaults;

import com.google.common.collect.ImmutableList;

/**
 * Greedy is good: An experimental study on minimum clique cover and maximum independent set problems for randomly generated rectangles
 * 
 * Finding the connected components and a maximum clique of an intersection graph of rectangles in the plane, Journal of Algorithms, Volume 4, Issue 4, December 1983, Pages 310–323
 * A note on maximum independent sets in rectangle intersection graphs, Information Processing Letters, Volume 89, Issue 1, 16 January 2004, Pages 19–23
 * 
 * Partitions a rectangle graph (boxicity=2) into cliques such that
 * each vertex is contained in exactly one clique.
 * 
 * Partitioning is performed greedily with nodes iteratively allocated to the
 * maximum clique.
 * 
 * This is similar to the classical maximum clique decomposition problem but
 * instead of removing maximum clique edges, maximum verticies are removed.
 *  
 * This algorithm is optimised for sparse rectangle graphs. 
 * 
 * @author Daniel Cameron
 */
public class RectangleGraphMaximalCliqueCalculator {
	private GraphNode lastNode = null;
	private List<GraphNode> outBuffer;
	private final PriorityQueue<GraphNode> activeEndingX = new PriorityQueue<GraphNode>(11, GraphNode.ByEndXStartYEndY); // sorted by endX
	private final PriorityQueue<GraphNode> activeScanlineStartingNodesEndingY = new PriorityQueue<GraphNode>(11, GraphNode.ByEndY); // sorted by endX
	private final ScanlineInterval activeScanlineStart = new ScanlineInterval();
	private ScanlineInterval activeScanlineCurrentPosition = activeScanlineStart;
	private int activeScanlineActiveWeight = 0;
	private long scanlineX = Long.MIN_VALUE;
	private class ScanlineInterval implements Cloneable {
		private ScanlineInterval(
				long startY,
				long startx,
				float weight,
				int startHere,
				int endHere,
				ScanlineInterval next) {
			this.startY = startY;
			this.startX = startx;
			this.weight = weight;
			this.startHere = startHere;
			this.endHere = endHere;
			this.next = next;
		}
		public ScanlineInterval() {
			this(Long.MIN_VALUE, Long.MIN_VALUE, 0, 0, 0, null);
		}
		/**
		 * Splits the given node at the given position
		 * @param y y value to split at
		 * @return node starting at the given y value
		 */
		public void splitAt(long y) {
			assert(y >= getStartY());
			assert(y <= getEndY());
			if (y == startY) {
				return;
			}
			ScanlineInterval newNode = new ScanlineInterval(
					y,
					this.startX,
					this.weight,
					0,
					this.endHere,
					this.next);
			this.endHere = 0;
			this.next = newNode;
		}
		private long startY = Long.MIN_VALUE;
		private long startX = Long.MIN_VALUE;
		private float weight = 0;
		private int startHere = 0;
		private int endHere = 0;
		private ScanlineInterval next;
		public long getStartY() {
			return startY;
		}
		public long getEndY() {
			if (next == null) return Long.MAX_VALUE;
			return next.startY - 1;
		}
		private boolean isMaximalClique() {
			return startHere > 0 && endHere > 0;
		}
		@Override
		public String toString() {
			return String.format("x=%d y=%d start#=%d end#=%d %f", startX, startY, startHere, endHere, weight);
		}
		public ScanlineInterval getNext() {
			return next;
		}
		public void mergeWithNext() { 
			assert(weight == next.weight);
			assert(next != null);
			next = next.next;
		}
	}
	private List<GraphNode> getCalledCliques() {
		List<GraphNode> result = outBuffer == null ? ImmutableList.<GraphNode>of() : outBuffer;
		outBuffer = null;
		return result;
	}
	/**
	 * Advances to the next position
	 * @param node
	 * @return
	 */
	public List<GraphNode> next(GraphNode node) {
		assert(lastNode == null || GraphNode.ByStartXY.compare(lastNode, node) <= 0);
		assert(node.startX >= scanlineX);
		lastNode = node;
		if (node.startX != scanlineX) {
			// finish off current scanline
			scanlineProcessYEndBefore(Long.MAX_VALUE, false);
			assert(activeScanlineStartingNodesEndingY.isEmpty());
			// advance scanline
			processEndXBefore(node.startX);
			// reset at start for this new scanline
			activeScanlineCurrentPosition = activeScanlineStart;
			scanlineX = node.startX;
		}
		incorporateInCurrentScanline(node, false);
		if (Defaults.PERFORM_EXPENSIVE_CLIQUE_SANITY_CHECKS) sanityCheck();
		return getCalledCliques();
	}
	private void mergeNodes() {
		assert(activeScanlineStartingNodesEndingY.isEmpty()); // can't merge in the middle of processing
		for (ScanlineInterval si = activeScanlineStart; si != null && si.next != null; si = si.next) {
			if (si.endHere == 0 && si.next.startHere == 0) {
				// we can merge these together since the separating node is no longer around
				si.mergeWithNext();
			}
		}
	}
	private void sanityCheck() {
		// check activeScanline is order
		for (ScanlineInterval lastsi = null, si = activeScanlineStart; si != null; lastsi = si, si = si.next) {
			assert(si.getStartY() <= si.getEndY());
			if (lastsi != null) {
				assert(lastsi.getStartY() < si.getStartY());
				assert(lastsi.getEndY() == si.getStartY() - 1);
			} else {
				assert(si.getStartY() == Long.MIN_VALUE);
			}
			if (si.next == null) {
				assert(si.getEndY() == Long.MAX_VALUE);
			}
		}
	}
	/**
	 * Advances the current scanline position to the starting interval
	 * of the given node and adds the given node to the scanline active
	 * set 
	 * @param node node to start processing
	 */
	private void incorporateInCurrentScanline(GraphNode node, boolean removingNodes) {
		assert(scanlineX == node.startX);
		assert(activeScanlineCurrentPosition.startY <= node.startY);
		scanlineProcessYEndBefore(node.startY, removingNodes);
		activeScanlineCurrentPosition.splitAt(node.startY);
		advanceScanlineToIntervalContaining(node.startY, removingNodes);
		if (removingNodes) {
			activeScanlineActiveWeight -= node.weight;
			activeScanlineCurrentPosition.startHere--;
		} else {
			activeScanlineActiveWeight += node.weight;
			activeScanlineCurrentPosition.startHere++;
			activeScanlineCurrentPosition.startX = node.startX;
		}
		activeScanlineStartingNodesEndingY.add(node);
		assert(activeScanlineStartingNodesEndingY.peek().endY >= node.startY);
		assert(activeScanlineCurrentPosition.startY == node.startY);
		if (Defaults.PERFORM_EXPENSIVE_CLIQUE_SANITY_CHECKS) assert(activeScanlineStartingNodesEndingY.contains(node));
	}
	/**
	 * Advances the current scanline position to the given y position.
	 * @param endYBefore position to advance to
	 * @param callCliques call maximal cliques
	 */
	private void scanlineProcessYEndBefore(long endYBefore, boolean removingNodes) {
		while (!activeScanlineStartingNodesEndingY.isEmpty() && activeScanlineStartingNodesEndingY.peek().endY < endYBefore) {
			GraphNode completed = activeScanlineStartingNodesEndingY.poll();
			activeScanlineActiveWeight -= completed.weight;
			advanceScanlineToIntervalContaining(completed.endY, removingNodes);
			if (activeScanlineCurrentPosition.getEndY() != completed.endY) { 
				activeScanlineCurrentPosition.splitAt(completed.endY + 1);
			}
			assert(activeScanlineCurrentPosition.getEndY() == completed.endY);
			if (removingNodes) {
				activeScanlineCurrentPosition.endHere--;
				activeScanlineCurrentPosition.weight -= activeScanlineActiveWeight;
			} else {
				activeScanlineCurrentPosition.endHere++;
				activeScanlineCurrentPosition.weight += activeScanlineActiveWeight;
			}
		}
		// advance position to node containing endYBefore
		advanceScanlineToIntervalContaining(endYBefore, removingNodes);
	}
	/**
	 * Advances the current scanline to the interval containing
	 * @param y
	 * @return 
	 */
	private ScanlineInterval advanceScanlineToIntervalContaining(long y, boolean removingNodes) {
		assert(activeScanlineCurrentPosition.getStartY() <= y); // can't advance backwards
		while (activeScanlineCurrentPosition.getEndY() < y) {
			if (removingNodes && activeScanlineCurrentPosition.isMaximalClique() && !activeScanlineStartingNodesEndingY.isEmpty()) {
				outBuffer.add(new GraphNode(
						activeScanlineCurrentPosition.startX, scanlineX,
						activeScanlineCurrentPosition.getStartY(), activeScanlineCurrentPosition.getEndY(),
						activeScanlineCurrentPosition.weight));
			}
			if (removingNodes) {
				activeScanlineCurrentPosition.weight -= activeScanlineActiveWeight;
			} else {
				activeScanlineCurrentPosition.weight += activeScanlineActiveWeight;
			}
			activeScanlineCurrentPosition = activeScanlineCurrentPosition.getNext();
		}
		assert(activeScanlineCurrentPosition.getStartY() <= y);
		assert(activeScanlineCurrentPosition.getEndY() >= y);
		if (y == Long.MAX_VALUE) {
			// if we hit our sentinel, we should be done processing
			assert(activeScanlineStartingNodesEndingY.isEmpty());
			assert(activeScanlineActiveWeight == 0);
			assert(activeScanlineCurrentPosition.next == null);
		}
		return activeScanlineCurrentPosition;
	}
	private void processEndXBefore(long endBeforeX) {
		outBuffer = new ArrayList<GraphNode>();
		while (!activeEndingX.isEmpty() && activeEndingX.peek().endX < endBeforeX) {
			scanlineX = activeEndingX.peek().endX;
			processEndingXOnCurrentScanline();
		}
	}
	private void processEndingXOnCurrentScanline() {
		activeScanlineCurrentPosition = activeScanlineStart;
		while (activeEndingX.peek().endX == scanlineX) {
			GraphNode completion = activeEndingX.poll();
			incorporateInCurrentScanline(completion, true);
		}
		scanlineProcessYEndBefore(Long.MAX_VALUE, true);
		mergeNodes();
	}
	public List<GraphNode> complete() {
		processEndXBefore(Long.MAX_VALUE);
		return outBuffer;
	}
}