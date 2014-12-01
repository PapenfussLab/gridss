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
	//private final RangeMap<GraphNode> activeScanlineEvidence;
	/**
	 * Contains GraphNodes of which the start Y has been processed but the end Y has not yet been encountered
	 */
	private final PriorityQueue<GraphNode> activeScanlineEndingY = new PriorityQueue<GraphNode>(11, GraphNode.ByEndY); // sorted by endX
	private final ScanlineInterval activeScanlineStart;
	private ScanlineInterval activeScanlineCurrentPosition = null;
	private int activeScanlineActiveWeight = 0;
	private long scanlineX = Long.MIN_VALUE;
	public RectangleGraphMaximalCliqueCalculator() {
		activeScanlineStart = new ScanlineInterval();
		activeScanlineStart.next = new ScanlineInterval();
		activeScanlineStart.next.startY = Long.MAX_VALUE;
		activeScanlineCurrentPosition = activeScanlineStart;
	}
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
		 * Splits the given node containing the given y value
		 * so a node starts at the given position
		 * @param y start y to ensure
		 */
		public void splitStartAt(long y) {
			assert(y >= getStartY());
			assert(y <= getEndY());
			if (y == getStartY()) return;
			ScanlineInterval newNode = new ScanlineInterval(
					y,
					this.startX,
					this.weight,
					0,
					this.endHere,
					this.next);
			this.endHere = 0;
			this.next = newNode;
			this.couldBeMaximalClique = false;
		}
		/**
		 * Splits the given node containing the given y value
		 * so a node ends at the given position
		 * @param y end y to ensure
		 */
		public void splitEndAt(long y) {
			if (getEndY() == y) return;
			splitStartAt(y + 1);
		}
		private long startY = Long.MIN_VALUE;
		private long startX = Long.MIN_VALUE;
		private float weight = 0;
		private int startHere = 0;
		private int endHere = 0;
		private ScanlineInterval next;
		private boolean couldBeMaximalClique = false;
		public long getStartY() {
			return startY;
		}
		public long getEndY() {
			if (next == null) return Long.MAX_VALUE;
			return next.startY - 1;
		}
		private boolean isMaximalClique() {
			return startHere > 0 && endHere > 0 && couldBeMaximalClique;
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
			scanlineCompleteProcessing(1);
			// advance scanline
			processEndXBefore(node.startX);
			scanlineX = node.startX;
		}
		incorporateInCurrentScanline(node, 1);
		activeEndingX.add(node);
		return getCalledCliques();
	}
	private void mergeIntervals() {
		assert(activeScanlineEndingY.isEmpty()); // can't merge in the middle of processing
		for (ScanlineInterval si = activeScanlineStart; si != null && si.next != null; si = si.next) {
			if (si.endHere == 0 && si.next.startHere == 0 && si.next.getStartY() != Long.MAX_VALUE) {
				// we can merge these together since the separating node is no longer around
				si.mergeWithNext();
			}
		}
	}
	private boolean sanityCheckScanlineActive() {
		if (!Defaults.PERFORM_EXPENSIVE_CLIQUE_SANITY_CHECKS) return true;
		assert(activeScanlineActiveWeight > 0);
		assert(!activeScanlineEndingY.isEmpty()); 
		assert(activeScanlineCurrentPosition != null);
		assert(activeScanlineCurrentPosition.getStartY() < Long.MAX_VALUE);
		assert(activeScanlineCurrentPosition.getEndY() < Long.MAX_VALUE);
		assert(sanityCheck());
		return true;
	}
	private boolean sanityCheckScanlineComplete() {
		if (!Defaults.PERFORM_EXPENSIVE_CLIQUE_SANITY_CHECKS) return true;
		assert(activeScanlineCurrentPosition == activeScanlineStart);
		assert(activeScanlineEndingY.isEmpty());
		assert(activeScanlineActiveWeight == 0);
		for (ScanlineInterval si = activeScanlineStart.next; si != null && si.next != null; si = si.next) {
			if (si.next.next != null) {
				// something should be splitting this node from the previous one
				assert(si.endHere > 0 || si.next.startHere > 0);
			}
		}
		assert(sanityCheck());
		return true;
	}
	private boolean sanityCheck() {
		assert(activeScanlineStart != null);
		assert(activeScanlineStart.getStartY() == Long.MIN_VALUE);
		assert(activeScanlineStart.weight == 0);
		assert(activeScanlineStart.startHere == 0);
		assert(activeScanlineStart.endHere == 0);
		// check activeScanline is ordered
		for (ScanlineInterval lastsi = null, si = activeScanlineStart; si != null; lastsi = si, si = si.next) {
			assert(si.getStartY() <= si.getEndY());
			if (lastsi != null) {
				assert(lastsi.getStartY() < si.getStartY());
				assert(lastsi.getEndY() == si.getStartY() - 1);
			} else {
				assert(si.getStartY() == Long.MIN_VALUE);
			}
			if (si.next == null) {
				// sentinel
				assert(si.getStartY() == Long.MAX_VALUE);
				assert(si.getEndY() == Long.MAX_VALUE);
				assert(si.weight == 0);
				assert(si.startHere == 0);
				assert(si.endHere == 0);
			}
		}
		float weight = 0;
		for (GraphNode n : activeScanlineEndingY) {
			assert(n.startY <= activeScanlineCurrentPosition.getStartY());
			weight += n.weight;
		}
		assert(activeScanlineActiveWeight == weight);
		return true;
	}
	/**
	 * Advances the current scanline position to the starting interval
	 * of the given node and adds the given node to the scanline active
	 * set 
	 * @param node node to start processing
	 * @Param multiplier 1 indicates we are incorporating the start of the given GraphNode to the current scanline 
	 * 0 indicates we are incorporating the end of the given GraphNode to the current scanline
	 */
	private void incorporateInCurrentScanline(GraphNode node, int multiplier) {
		assert(multiplier == -1 || multiplier == 1);
		assert(scanlineX == (multiplier == 1 ? node.startX : node.endX));
		long y = node.startY;
		assert(activeScanlineCurrentPosition.startY <= y);
		scanlineProcessYEndBefore(y, multiplier);
		activeScanlineCurrentPosition.splitStartAt(y);
		scanlineProcessYEndBefore(y, multiplier);
		activeScanlineActiveWeight += node.weight;
		activeScanlineCurrentPosition.startHere += multiplier;
		if (multiplier == 1) {
			activeScanlineCurrentPosition.startX = node.startX;
		}
		activeScanlineEndingY.add(node);
		assert(activeScanlineCurrentPosition.getStartY() == y);
		assert(activeScanlineEndingY.contains(node));
		assert(sanityCheckScanlineActive());
	}
	/**
	 * Advances the current scanline position to the given y position.
	 * @param endYBefore position to advance to
	 * @Param multiplier 1 indicates we are incorporating the start of the given GraphNode to the current scanline 
	 * 0 indicates we are incorporating the end of the given GraphNode to the current scanline
	 */
	private void scanlineProcessYEndBefore(long endYBefore, int multiplier) {
		assert(multiplier == -1 || multiplier == 1);
		while (!activeScanlineEndingY.isEmpty() && activeScanlineEndingY.peek().endY < endYBefore) {
			GraphNode node = activeScanlineEndingY.poll();
			long y = node.endY;
			int yendCount = 1;
			float yendWeight = node.weight;
			while (!activeScanlineEndingY.isEmpty() && activeScanlineEndingY.peek().endY == y) {
				node = activeScanlineEndingY.poll();
				yendCount++;
				yendWeight += node.weight;
			}
			advanceScanlineToIntervalContaining(y, multiplier);
			activeScanlineCurrentPosition.splitEndAt(y);
			activeScanlineCurrentPosition.endHere += yendCount * multiplier;
			advanceScanlineToIntervalContaining(y + 1, multiplier);
			activeScanlineActiveWeight -= yendWeight;
		}
		// advance position to node containing endYBefore
		advanceScanlineToIntervalContaining(endYBefore, multiplier);
	}
	/**
	 * Advances the current scanline to the interval containing
	 * This method is responsible for setting activeScanlineCurrentPosition
	 * This method is responsible for updating scanline weights based on activeScanlineActiveWeight
	 * @param y
	 * @return 
	 */
	private void advanceScanlineToIntervalContaining(long y, int multiplier) {
		assert(activeScanlineCurrentPosition.getStartY() <= y); // can't advance backwards
		while (activeScanlineCurrentPosition.getEndY() < y) {
			activeScanlineCurrentPosition.weight += activeScanlineActiveWeight * multiplier;
			if (activeScanlineActiveWeight != 0) {
				// could be maximal if we're adding new evidence
				// if we're removing evidence then we're now definitely not maximal
				// if we're doing neither then there is no change from the previous scanline
				activeScanlineCurrentPosition.couldBeMaximalClique = multiplier == 1;
			}
			activeScanlineCurrentPosition = activeScanlineCurrentPosition.getNext();
		}
		assert(activeScanlineCurrentPosition.getStartY() <= y);
		assert(activeScanlineCurrentPosition.getEndY() >= y);
	}
	/**
	 * Calls maximum cliques
	 * @param endingCurrentScanline nodes ending here. Maximum cliques will always occur within one of these intervals
	 */
	private void callMaximumCliques(List<GraphNode> endingCurrentScanline) {
		ScanlineInterval interval = activeScanlineStart;
		int index = 0;
		while (index < endingCurrentScanline.size()) {
			long startY = endingCurrentScanline.get(index).startY;
			long endY = endingCurrentScanline.get(index).endY;
			index++;
			while (index + 1 < endingCurrentScanline.size() && endingCurrentScanline.get(index + 1).startY <= endY) {
				// expand the current calling interval due to overlap
				endY = Math.max(endY, endingCurrentScanline.get(index + 1).endY);
				index++;
			}
			// advance to interval
			while (interval.getEndY() < startY) {
				assert(interval.next != null);
				interval = interval.next;
			}
			// call cliques in interval
			assert(interval.getStartY() == startY);
			while (interval.getStartY() <= endY) {
				if (interval.isMaximalClique()) {
					outBuffer.add(new GraphNode(
							interval.startX, scanlineX,
							interval.getStartY(), interval.getEndY(),
							interval.weight));
				}
				interval = interval.next;
			}
			assert(interval.getStartY() == endY + 1);
		}
	}
	private void scanlineCompleteProcessing(int multiplier) {
		scanlineProcessYEndBefore(Long.MAX_VALUE, multiplier);
		// reset ready for next scanline
		activeScanlineCurrentPosition = activeScanlineStart;
		if (multiplier == -1) {
			// removal of nodes can result in adjacent intervals requiring merge
			mergeIntervals();
		}
		assert(sanityCheckScanlineComplete());
	}
	private void processEndXBefore(long endBeforeX) {
		outBuffer = new ArrayList<GraphNode>();
		while (!activeEndingX.isEmpty() && activeEndingX.peek().endX < endBeforeX) {
			scanlineX = activeEndingX.peek().endX;
			processEndingXOnCurrentScanline();
		}
	}
	private void processEndingXOnCurrentScanline() {
		assert(activeScanlineEndingY.isEmpty());
		List<GraphNode> endingCurrentScanline = new ArrayList<GraphNode>();
		while (!activeEndingX.isEmpty() && activeEndingX.peek().endX == scanlineX) {
			endingCurrentScanline.add(activeEndingX.poll());
		}
		callMaximumCliques(endingCurrentScanline);
		for (GraphNode g : endingCurrentScanline) {
			incorporateInCurrentScanline(g, -1);
		}
		scanlineCompleteProcessing(-1);
	}
	public List<GraphNode> complete() {
		scanlineCompleteProcessing(1);
		processEndXBefore(Long.MAX_VALUE);
		return outBuffer;
	}
}