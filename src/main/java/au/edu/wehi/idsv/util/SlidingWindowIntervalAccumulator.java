package au.edu.wehi.idsv.util;

import java.util.NavigableMap;
import java.util.PriorityQueue;
import java.util.TreeMap;

/**
 * Counts the number of intervals added for a given position.
 * 
 * Intervals are counted inclusive of both end positions [start, end]
 * 
 * For any given position, all intervals overlapping or preceding the
 * given position must be added to the accumulator before invoking
 * getIntervalCount() on that position.
 * 
 * Results for positions before a given sliding window are discarded. 
 * 
 * @author Daniel Cameron
 *
 */
public class SlidingWindowIntervalAccumulator implements Cloneable {
	/**
	 * Get the level of reference support (reads and concordant pairs) for the lack
	 * of breakpoint immediately after the given position 
	 * 
	 * @param position
	 * @return number of reads & concordant pairs supporting reference
	 */
	public int getIntervalCount(int position) {
		if (position > currentPosition) {
			calculateTo(position);
		}
		Integer count = resultCache.get(position);
		if (count == null) throw new IllegalArgumentException(String.format("Query position (%d) not available. Query position must be (%d) or higher", position, currentPosition));
		return count;
	}
	// Optimisation: advance to min(startInterval.peek(), endInterval.peek()) with current support
	private void calculateTo(int position) {
		while (currentPosition < position) {
			currentPosition++;
			// start adds support
			while (!startInterval.isEmpty() && startInterval.peek() <= currentPosition) {
				startInterval.poll();
				currentSupport++;
			}
			// end removes support
			while (!endInterval.isEmpty() && endInterval.peek() < currentPosition) {
				endInterval.poll();
				currentSupport--;
			}
			resultCache.put(currentPosition, currentSupport);
		}
	}
	/**
	 * Explicitly advances the window of validity by discarding all interval counts
	 * preceding the given position
	 * @param earliest position to retain a valid result for
	 */
	public void advanceTo(int position) {
		// Advance calculation if necessary, discarding results as they are outside
		// the window of validity
		if (currentPosition < position) {
			currentPosition = position;
			while (!startInterval.isEmpty() && startInterval.peek() <= currentPosition) {
				startInterval.poll();
				currentSupport++;
			}
			while (!endInterval.isEmpty() && endInterval.peek() < currentPosition) {
				endInterval.poll();
				currentSupport--;
			}
		}
		while (!resultCache.isEmpty() && resultCache.firstKey() < position) {
			resultCache.remove(resultCache.firstKey());
		}
	}
	/**
	 * Adds a new interval to the accumulator.
	 * 
	 * @param start inclusive start of interval
	 * @param end inclusive end of interval
	 */
	public void addInterval(int start, int end) {
		if (end < start) throw new IllegalArgumentException(String.format("end of interval (%d) cannot preceed start (%d)", end, start));
		if (start <= currentPosition) throw new IllegalArgumentException(String.format("Cannot add interval starting at %d when accumulator has processed up to %d", start, currentPosition));
		startInterval.add(start);
		endInterval.add(end);
	}
	// Optimisation: fast advancement when 0 support/interval list(s) are empty
	private int currentPosition = -1;
	private int currentSupport = 0;
	private PriorityQueue<Integer> startInterval = new PriorityQueue<Integer>();
	private PriorityQueue<Integer> endInterval = new PriorityQueue<Integer>();
	// Optimisation: find/write more efficient data structure (CircularArrayList or IntervalTree)
	private NavigableMap<Integer, Integer> resultCache = new TreeMap<Integer, Integer>();
	@Override
	public SlidingWindowIntervalAccumulator clone() {
		SlidingWindowIntervalAccumulator result = new SlidingWindowIntervalAccumulator();
		result.currentPosition = this.currentPosition;
		result.currentSupport = this.currentSupport;
		result.startInterval = new PriorityQueue<Integer>(this.startInterval);
		result.endInterval = new PriorityQueue<Integer>(this.endInterval);
		result.resultCache = new TreeMap<Integer, Integer>(this.resultCache);
		return result;
	}
}
