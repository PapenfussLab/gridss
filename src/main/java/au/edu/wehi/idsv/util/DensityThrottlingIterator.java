package au.edu.wehi.idsv.util;

import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.PeekingIterator;
import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import it.unimi.dsi.fastutil.longs.LongPriorityQueue;

import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * Randomly filters sorted records to limit the record density
 * 
 * All records are accepted up to a threshold value, then
 * a exponential back-off is used to randomly filter additional
 * records to ensure an average maximum number of records in
 * any given window. 
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public abstract class DensityThrottlingIterator<T> implements PeekingIterator<T>, TrackedBuffer {
	private final Iterator<T> underlying;
	private final double windowSize;
	private final double acceptDensity;
	private final double maxDensity;
	private final LongPriorityQueue inWindow = new LongArrayFIFOQueue();
	private final LongPriorityQueue emittedInWindow = new LongArrayFIFOQueue();
	private final Random random = new Random(0); // Seed set for reproducible results
	private T nextRecord = null;
	
	/**
	 * @param it iterator to filter. Cannot contain null elements
	 * @param windowSize Size of window to track density over
	 * @param acceptAllCount number of record in window before throttling starts
	 * @param maxCount maximum average number of records in window 
	 */
	public DensityThrottlingIterator(Iterator<T> it, int windowSize, double acceptDensity, double targetDensity) {
		this.underlying = it;
		this.windowSize = windowSize;
		this.acceptDensity = acceptDensity;
		this.maxDensity = targetDensity;
	}
	protected abstract long getPosition(T record);
	protected abstract boolean excludedFromThrottling(T record);
	private void ensureNext() {
		while (nextRecord == null && underlying.hasNext()) {
			nextRecord = underlying.next();
			if (excludedFromThrottling(nextRecord)) {
				return;
			}
			long position = getPosition(nextRecord);
			// remove records
			while (!inWindow.isEmpty() && inWindow.firstLong() <= position - windowSize) {
				inWindow.dequeueLong();
			}
			while (!emittedInWindow.isEmpty() && emittedInWindow.firstLong() <= position - windowSize) {
				emittedInWindow.dequeueLong();
			}
			inWindow.enqueue(position);
			if (!isFiltered(position, nextRecord)) {
				emittedInWindow.enqueue(position);
			} else {
				nextRecord = null;
			}
		}
	}
	private boolean isFiltered(long position, T record) {
		if (isBelowUnconditionalAcceptanceThreshold()) {
			// accept all record under the threshold
			return false;
		}
		double x = ((inWindow.size() / windowSize) - acceptDensity) / maxDensity;
		if (Math.exp(-x) >= random.nextDouble()) {
			// exponential back-off did not filter
			emittedInWindow.enqueue(position);
			return false;
		}
		return true;
	}
	public boolean isBelowUnconditionalAcceptanceThreshold() {
		return emittedInWindow.size() / windowSize < acceptDensity;
	}
	public double currentDensity() {
		return inWindow.size() / windowSize;
	}
	@Override
	public boolean hasNext() {
		ensureNext();
		return nextRecord != null;
	}
	@Override
	public T next() {
		ensureNext();
		T r = nextRecord;
		nextRecord = null;
		return r;
	}
	@Override
	public T peek() {
		ensureNext();
		return nextRecord;
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	private String trackedBufferName_inWindow = "densityThrottle.inWindow";
	private String trackedBufferName_emittedInWindow = "densityThrottle.emittedInWindow";
	@Override
	public void setTrackedBufferContext(String context) {
		this.trackedBufferName_inWindow = context + ".densityThrottle.inWindow";
		this.trackedBufferName_emittedInWindow = context + ".densityThrottle.emittedInWindow";
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return ImmutableList.of(
				new NamedTrackedBuffer(trackedBufferName_inWindow, inWindow.size()),
				new NamedTrackedBuffer(trackedBufferName_emittedInWindow, emittedInWindow.size())
				);
	}
}
