package au.edu.wehi.idsv.util;

import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import com.google.common.base.Function;
import com.google.common.collect.*;
import com.google.common.primitives.Longs;
import htsjdk.samtools.util.Log;

import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Sorts a mostly-sorted input sequence.
 * 
 * As SAM/BAM input is sorted by alignment start position, sorting on evidence
 * position does not require a full sort as the difference between breakend
 * start position and the alignment start position is bounded by the fragment size
 * for read pair evidence, and the read length for soft clip evidence.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class WindowedSortingIterator<T> extends AbstractIterator<T> implements TrackedBuffer {
	private static final Log log = Log.getInstance(WindowedSortingIterator.class);
	private final PriorityQueue<T> calls; 
	private final long windowSize;
	private final PeekingIterator<T> it;
	private final Function<T, Long> toCoordinate;
	private long lastPosition = Long.MIN_VALUE;
	private final Comparator<T> sortOrder;
	private T lastEmitted = null;
	/**
	 * Creates a new sorted iterator from a mostly-sorted sequence
	 * @param it mostly-sorted sequence. Records cannot be out of order by more than windowSize
	 * @param transform Coordinate transform for position of record.
	 * @param windowSize Maximum coordinate-space length that records can deviate from a sorted sequence 
	 */
	public WindowedSortingIterator(final Iterator<T> it, final Function<T, Long> transform, final long windowSize) {
		this(it, transform, windowSize, new Ordering<T>() {
			public int compare(T arg0, T arg1) {
				return Longs.compare(transform.apply(arg0), transform.apply(arg1));
			  }
		});
	}
	public WindowedSortingIterator(final Iterator<T> it, final Function<T, Long> transform, final long windowSize, final Comparator<T> sortOrder) {
		this.windowSize = windowSize;
		this.it = Iterators.peekingIterator(it);
		this.toCoordinate = transform;
		this.calls = new PriorityQueue<T>(32, sortOrder);
		this.sortOrder = sortOrder;
	}
	@Override
	protected T computeNext() {
		advanceUnderlying();
		if (calls.isEmpty()) return endOfData();
		T next = calls.poll();
		long nextPos = toCoordinate.apply(next);
		if (nextPos < lastPosition) {
			log.error("Sanity check failure: sorting window size too small: evidence out of order at linear coordinate" + nextPos);
		}
		if (lastEmitted != null && sortOrder.compare(lastEmitted, next) > 0) {
			throw new IllegalStateException(String.format("Unable to sort output with window size of %d. %s emitted before %s", windowSize, lastEmitted, next));
		}
		lastEmitted = next;
		return next;
	}
	private void advanceUnderlying() {
		while (it.hasNext() && (calls.isEmpty() || nextRecordCouldBeAtStartOfWindow())) {
			T next = it.next();
			calls.add(next);
		}
	}
	private boolean nextRecordCouldBeAtStartOfWindow() {
		long bufferPosition = toCoordinate.apply(calls.peek());
		long nextPosition = toCoordinate.apply(it.peek());
		return nextPosition <= bufferPosition + windowSize;
	}
	private String trackedBufferName_calls = "windowedSort";
	@Override
	public void setTrackedBufferContext(String context) {
		this.trackedBufferName_calls = context + ".windowedSort";
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return ImmutableList.of(
				new NamedTrackedBuffer(trackedBufferName_calls, calls.size())
				);
	}
}
