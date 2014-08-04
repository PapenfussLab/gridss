package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.base.Function;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;
import com.google.common.primitives.Longs;

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
public abstract class WindowedSortingIterator<T> extends AbstractIterator<T> {
	private static final Log log = Log.getInstance(WindowedSortingIterator.class);
	private final PriorityQueue<T> calls; 
	private final long windowSize;
	private final PeekingIterator<T> it;
	private final Function<T, Long> toCoordinate;
	private long lastPosition = Long.MIN_VALUE;
	/**
	 * Creates a new sorted iterator from a mostly-sorted sequence
	 * @param it mostly-sorted sequence. Records cannot be out of order by more than windowSize
	 * @param transform Coordinate transform for position of record.
	 * @param windowSize Maximum coordinate-space length that records can deviate from a sorted sequence 
	 */
	public WindowedSortingIterator(Iterator<T> it, Function<T, Long> transform, long windowSize) {
		this.windowSize = windowSize;
		this.it = Iterators.peekingIterator(it);
		this.toCoordinate = transform;
		this.calls = new PriorityQueue<T>(32, new Ordering<T>() {
			public int compare(T arg0, T arg1) {
				return Longs.compare(toCoordinate.apply(arg0), toCoordinate.apply(arg1));
			  }
		});
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
		return next;
	}
	private void advanceUnderlying() {
		while (it.hasNext() && (calls.isEmpty() || nextRecordCouldBeAtStartOfWindow())) {
			T next = it.next();
			if (next == null) {
				throw new RuntimeException("Sanity check failure: null evidence");
			}
			calls.add(next);
		}
	}
	private boolean nextRecordCouldBeAtStartOfWindow() {
		long bufferPosition = toCoordinate.apply(calls.peek());
		long nextPosition = toCoordinate.apply(it.peek());
		return nextPosition <= bufferPosition + windowSize;
	}
}
