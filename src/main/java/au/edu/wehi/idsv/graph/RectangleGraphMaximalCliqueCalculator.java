package au.edu.wehi.idsv.graph;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.visualisation.TrackedState;
import com.google.common.collect.ImmutableList;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Calculates all maximal cliques of a rectangle graph
 *
 * Maximum clique problem of rectangle graphs, Advances in Computer Research, D T Lee, 1 (1983), pp. 91-107
 * 
 * Better implementation would be to partition the rectangle graph (boxicity=2)
 * into cliques such that each vertex is contained in exactly one clique.
 * 
 * Since total vertex weight is fixed, minimising number of cliques is equivalent to maximising average clique weight.
 * This an NP, see:
 * Greedy is good: An experimental study on minimum clique cover and maximum independent set problems for randomly generated rectangles
 * Finding the connected components and a maximum clique of an intersection graph of rectangles in the plane, Journal of Algorithms, Volume 4, Issue 4, December 1983, Pages 310-323
 * A note on maximum independent sets in rectangle intersection graphs, Information Processing Letters, Volume 89, Issue 1, 16 January 2004, Pages 19-23
 * GREEDY MAXIMUM-CLIQUE DECOMPOSITIONS http://faculty.tru.ca/smcguinness/greedymaxclique.pdf (we want to decompose by removing vertices, not edges) 
 * 
 * @author Daniel Cameron
 */
public class RectangleGraphMaximalCliqueCalculator implements TrackedState {
	private RectangleGraphNode lastNode = null;
	private List<RectangleGraphNode> outBuffer;
	private final PriorityQueue<RectangleGraphNode> activeEndingX = new PriorityQueue<RectangleGraphNode>(11, RectangleGraphNode.ByEndXStartYEndY); // sorted by endX
	//private final RangeMap<GraphNode> activeScanlineEvidence;
	/**
	 * Contains GraphNodes of which the start Y has been processed but the end Y has not yet been encountered
	 */
	private final PriorityQueue<RectangleGraphNode> activeScanlineEndingY = new PriorityQueue<RectangleGraphNode>(11, RectangleGraphNode.ByEndY); // sorted by endX
	private final ScanlineInterval activeScanlineStart;
	private ScanlineInterval activeScanlineCurrentPosition = null;
	private long activeScanlineActiveWeight = 0;
	private long scanlineX = Long.MIN_VALUE;
	public RectangleGraphMaximalCliqueCalculator() {
		activeScanlineStart = new ScanlineInterval();
		activeScanlineStart.next = new ScanlineInterval();
		activeScanlineStart.next.startY = Long.MAX_VALUE - 1;
		activeScanlineCurrentPosition = activeScanlineStart;
		assert(sanityCheckScanlineComplete());
	}

	/**
	 * Scanline interval of the rectangle graph.
	 * Scanline coordinates use half-open intervals.
	 * This differs from GraphNode representation
	 * @author Daniel Cameron
	 *
	 */
	private class ScanlineInterval implements Cloneable {
		private ScanlineInterval(
				long startY,
				long startx,
				long weight,
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
			this(Long.MIN_VALUE, Long.MAX_VALUE, 0, 0, 0, null);
		}
		/**
		 * Splits the given node containing the given y value
		 * so a node starts at the given position
		 * @param y start y to ensure
		 */
		public void splitAt(long y) {
			assert(y >= getStartY());
			assert(y < getEndY());
			if (y == getStartY()) return;
			ScanlineInterval newNode = new ScanlineInterval(
					y,
					Long.MAX_VALUE,
					this.weight,
					0,
					this.endHere,
					this.next);
			this.endHere = 0;
			this.next = newNode;
			this.startX = Long.MAX_VALUE;
		}
		private long startY = Long.MIN_VALUE;
		/**
		 * Long.MIN_VALUE indicates this interval is not maximal
		 */
		private long startX = Long.MAX_VALUE;
		private long weight = 0;
		private int startHere = 0;
		private int endHere = 0;
		private ScanlineInterval next;
		/**
		 * Start coordinate of the half-open interval
		 * @return
		 */
		public long getStartY() {
			return startY;
		}
		/**
		 * End coordinate of the half-open interval
		 * @return
		 */
		public long getEndY() {
			if (next == null) return Long.MAX_VALUE;
			return next.startY;
		}
		private boolean isMaximalClique() {
			return startX != Long.MAX_VALUE;
		}
		public ScanlineInterval getNext() {
			return next;
		}
		public void mergeWithNext() {
			assert(next != null);
			assert(next.next != null); // can't merge with the end sentinal
			assert(weight == next.weight);
			assert(startX == Long.MAX_VALUE);
			assert(next.startX == Long.MAX_VALUE);
			endHere = next.endHere;
			next = next.next;
		}
		@Override
		public String toString() {
			String s = String.format("[%d,%d)(w=%d,s=%d,e=%d,x=%d)", getStartY(), getEndY(), weight, startHere, endHere, startX);
			if (next != null) return s + "\n" + next.toString();
			return s;
		}
	}
	private List<RectangleGraphNode> getCalledCliques() {
		List<RectangleGraphNode> result = outBuffer == null ? ImmutableList.<RectangleGraphNode>of() : outBuffer;
		outBuffer = null;
		return result;
	}
	/**
	 * Advances to the next position
	 * @param node
	 * @return
	 */
	public List<RectangleGraphNode> next(RectangleGraphNode node) {
		assert(node.startX <= node.endX);
		assert(node.startY <= node.endY);
		assert(node.weight > 0);
		assert(node.startX >= scanlineX);
		assert(lastNode == null || RectangleGraphNode.ByStartXY.compare(lastNode, node) <= 0);
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
		// (make sure we don't merge our end sentinal
		for (ScanlineInterval si = activeScanlineStart; si != null && si.next != null && si.next.next != null; si = si.next) {
			// merge consecutive nodes
			while (si.endHere == 0 && si.next.startHere == 0 && si.next.next != null) {
				assert(si.weight == si.next.weight);
				// we can merge these together since the separating node is no longer around
				si.mergeWithNext();
			}
		}
	}
	private boolean sanityCheckScanlineActive() {
		if (!Defaults.SANITY_CHECK_CLIQUE) return true;
		assert(sanityCheck());
		assert(!activeScanlineEndingY.isEmpty());
		assert(activeScanlineActiveWeight > 0);
		assert(activeScanlineCurrentPosition != null);
		assert(activeScanlineCurrentPosition.getStartY() < Long.MAX_VALUE - 1);
		assert(activeScanlineCurrentPosition.getEndY() < Long.MAX_VALUE);
		return true;
	}
	private boolean sanityCheckScanlineComplete() {
		if (!Defaults.SANITY_CHECK_CLIQUE) return true;
		assert(sanityCheck());
		assert(activeScanlineCurrentPosition == activeScanlineStart);
		assert(activeScanlineEndingY.isEmpty());
		assert(activeScanlineActiveWeight == 0);
		for (ScanlineInterval si = activeScanlineStart.next; si != null && si.next != null; si = si.next) {
			if (si.next.next != null) {
				// something should be splitting this node from the previous one
				assert(si.endHere > 0 || si.next.startHere > 0);
			}
		}
		return true;
	}
	private boolean sanityCheck() {
		if (!Defaults.SANITY_CHECK_CLIQUE) return true;
		assert(activeScanlineStart != null);
		assert(activeScanlineStart.getStartY() == Long.MIN_VALUE);
		assert(activeScanlineStart.weight == 0);
		assert(activeScanlineStart.startHere == 0);
		assert(activeScanlineStart.endHere == 0);
		// check activeScanline is ordered
		for (ScanlineInterval lastsi = null, si = activeScanlineStart; si != null; lastsi = si, si = si.next) {
			assert(si.getStartY() < si.getEndY());
			if (lastsi != null) {
				assert(lastsi.getStartY() < si.getStartY());
				assert(lastsi.getEndY() == si.getStartY());
			} else {
				assert(si.getStartY() == Long.MIN_VALUE);
			}
			if (si.next == null) {
				// sentinel
				assert(si.getStartY() == Long.MAX_VALUE - 1);
				assert(si.getEndY() == Long.MAX_VALUE);
				assert(si.weight == 0);
				assert(si.startHere == 0);
				assert(si.endHere == 0);
			}
		}
		long weight = 0;
		for (RectangleGraphNode n : activeScanlineEndingY) {
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
	private void incorporateInCurrentScanline(RectangleGraphNode node, int multiplier) {
		assert(multiplier == -1 || multiplier == 1);
		assert(scanlineX == (multiplier == 1 ? node.startX : node.endX));
		long y = node.startY;
		assert(activeScanlineCurrentPosition.startY <= y);
		scanlineProcessYEndBefore(y, multiplier);
		if (activeScanlineCurrentPosition.getStartY() != y) {
			activeScanlineCurrentPosition.splitAt(y);
			scanlineProcessYEndBefore(y, multiplier);
		}
		activeScanlineActiveWeight += node.weight;
		activeScanlineCurrentPosition.startHere += multiplier;
		activeScanlineEndingY.add(node);
		assert(activeScanlineCurrentPosition.getStartY() == y);
		assert(activeScanlineEndingY.contains(node));
		assert(sanityCheckScanlineActive());
	}
	/**
	 * Advances the current scanline position to the half-open interval containing the given y position.
	 * @param endYBefore position to advance to
	 * @Param multiplier 1 indicates we are incorporating the start of the given GraphNode to the current scanline 
	 * 0 indicates we are incorporating the end of the given GraphNode to the current scanline
	 */
	private void scanlineProcessYEndBefore(long endYBefore, int multiplier) {
		assert(multiplier == -1 || multiplier == 1);
		while (!activeScanlineEndingY.isEmpty() && activeScanlineEndingY.peek().endY < endYBefore) {
			RectangleGraphNode node = activeScanlineEndingY.poll();
			long endYexclusive = node.endY + 1;
			int yendCount = 1;
			long yendWeight = node.weight;
			while (!activeScanlineEndingY.isEmpty() && activeScanlineEndingY.peek().endY + 1 == endYexclusive) {
				node = activeScanlineEndingY.poll();
				yendCount++;
				yendWeight += node.weight;
			}
			advanceScanlineToIntervalContaining(endYexclusive - 1, multiplier);
			if (activeScanlineCurrentPosition.getEndY() > endYexclusive) {
				activeScanlineCurrentPosition.splitAt(endYexclusive);
			}
			// no need to advance here since our current position is correct
			activeScanlineCurrentPosition.endHere += yendCount * multiplier;
			advanceScanlineToIntervalContaining(endYexclusive, multiplier); // move on past our closing position
			activeScanlineActiveWeight -= yendWeight;
		}
		// advance position to node containing endYBefore
		advanceScanlineToIntervalContaining(endYBefore, multiplier);
	}
	/**
	 * Advances the current scanline to the interval containing
	 * This method is responsible for setting activeScanlineCurrentPosition
	 * This method is responsible for updating scanline weights based on activeScanlineActiveWeight
	 * @param y included in half-open interval to advance scanline to  
	 */
	private void advanceScanlineToIntervalContaining(long y, int multiplier) {
		assert(activeScanlineCurrentPosition.getStartY() <= y); // can't advance backwards
		while (activeScanlineCurrentPosition.getEndY() <= y) {
			activeScanlineCurrentPosition.weight += activeScanlineActiveWeight * multiplier;
			if (activeScanlineActiveWeight != 0) {
				// could be maximal if we're adding new evidence
				// if we're removing evidence then we're now definitely not maximal
				// if we're doing neither then there is no change from the previous scanline
				activeScanlineCurrentPosition.startX = Long.MAX_VALUE;
				if (multiplier == 1 && activeScanlineCurrentPosition.startHere > 0 && activeScanlineCurrentPosition.endHere > 0) {
					activeScanlineCurrentPosition.startX = scanlineX;
				}
			}
			activeScanlineCurrentPosition = activeScanlineCurrentPosition.getNext();
		}
		assert(activeScanlineCurrentPosition.getStartY() <= y);
		assert(activeScanlineCurrentPosition.getEndY() > y);
	}
	/**
	 * Calls maximum cliques
	 * @param endingCurrentScanline nodes ending here. Maximum cliques will always occur within one of these intervals
	 */
	private void callMaximumCliques(List<RectangleGraphNode> endingCurrentScanline) {
		ScanlineInterval interval = activeScanlineStart;
		int index = 0;
		while (index < endingCurrentScanline.size()) {
			long startY = endingCurrentScanline.get(index).startY;
			long endYexclusive = endingCurrentScanline.get(index).endY + 1;
			index++;
			while (index < endingCurrentScanline.size() && endingCurrentScanline.get(index).startY <= endYexclusive) {
				// expand the current calling interval due to overlap
				endYexclusive = Math.max(endYexclusive, endingCurrentScanline.get(index).endY + 1);
				index++;
			}
			// advance to interval
			while (interval.getEndY() <= startY) {
				assert(interval.next != null);
				interval = interval.next;
			}
			// call cliques in interval
			assert(interval.getStartY() == startY);
			while (interval.getStartY() < endYexclusive) {
				if (interval.isMaximalClique()) {
					outBuffer.add(new RectangleGraphNode(
							interval.startX, scanlineX,
							interval.getStartY(), interval.getEndY() - 1, // convert back from half-open to close interval
							interval.weight));
				}
				interval = interval.next;
			}
			assert(interval.getStartY() == endYexclusive);
		}
	}
	private void scanlineCompleteProcessing(int multiplier) {
		scanlineProcessYEndBefore(Long.MAX_VALUE - 1, multiplier);
		// reset ready for next scanline
		activeScanlineCurrentPosition = activeScanlineStart;
		if (multiplier == -1) {
			// removal of nodes can result in adjacent intervals requiring merge
			mergeIntervals();
		}
		assert(sanityCheckScanlineComplete());
	}
	private void processEndXBefore(long endBeforeX) {
		outBuffer = new ArrayList<RectangleGraphNode>();
		while (!activeEndingX.isEmpty() && activeEndingX.peek().endX < endBeforeX) {
			scanlineX = activeEndingX.peek().endX;
			processEndingXOnCurrentScanline();
		}
	}
	private void processEndingXOnCurrentScanline() {
		assert(activeScanlineEndingY.isEmpty());
		List<RectangleGraphNode> endingCurrentScanline = new ArrayList<RectangleGraphNode>();
		while (!activeEndingX.isEmpty() && activeEndingX.peek().endX == scanlineX) {
			endingCurrentScanline.add(activeEndingX.poll());
		}
		callMaximumCliques(endingCurrentScanline);
		for (RectangleGraphNode g : endingCurrentScanline) {
			incorporateInCurrentScanline(g, -1);
		}
		scanlineCompleteProcessing(-1);
	}
	public List<RectangleGraphNode> complete() {
		scanlineCompleteProcessing(1);
		processEndXBefore(Long.MAX_VALUE);
		return outBuffer;
	}

	@Override
	public String[] trackedNames() {
		return new String[] {
			"outBufferSize",
			"activeEndingXSize",
			"activeScanlineEndingYSize",
		};
	}

	@Override
	public Object[] trackedState() {
		return new Object[] {
				outBuffer == null ? 0 : outBuffer.size(),
				activeEndingX == null ? 0 : activeEndingX.size(),
				activeScanlineEndingY == null ? 0 : activeScanlineEndingY.size(),
		};
	}

	@Override
	public Collection<TrackedState> trackedObjects() {
		return ImmutableList.of(this);
	}
}