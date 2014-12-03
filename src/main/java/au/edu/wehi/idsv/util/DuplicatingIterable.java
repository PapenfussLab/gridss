package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import com.google.common.collect.AbstractIterator;

/**
 * Duplicates the given iterator.
 * 
 * This wrapper is thread-safe.
 * 
 * <b>Separate producer and consumer threads are required as
 * iterator calls block the calling thread when sufficiently
 * further ahead of other iterators.
 * </b>
 * @author Daniel Cameron
 *
 */
public class DuplicatingIterable<T> implements Iterable<T> {
	private static final Log log = Log.getInstance(DuplicatingIterable.class);
	private final Iterator<T> it;
	private final int bufferSize;
	private final List<DuplicatingIterableIterator> iterators = new ArrayList<DuplicatingIterableIterator>();
	/**
	 * Duplicates an iterator
	 * @param it underlying iterator
	 * @param maxIteratorDifference maximum number of records an iterator can traverse before blocking to wait
	 * for other iterators to catch up
	 */
	public DuplicatingIterable(Iterator<T> it, int maxIteratorDifference) {
		if (it == null) throw new IllegalArgumentException();
		if (maxIteratorDifference <= 0) throw new IllegalArgumentException("buffer size must be greater than zero.");
		this.it = it;
		this.bufferSize = maxIteratorDifference;
	}
	/**
	 * Creates a new iterator from the current position
	 */
	@Override
	public synchronized Iterator<T> iterator() {
		DuplicatingIterableIterator consumer = new DuplicatingIterableIterator();
		iterators.add(consumer);
		return consumer;
	}
	private synchronized void consumeNext(DuplicatingIterableIterator invokingConsumer) {
		if (!invokingConsumer.queue.isEmpty()) {
			// another thread populated our queue while we were locked out attempting to do it ourself
			return;
		}
		if (it.hasNext()) {
			T n = it.next();
			for (DuplicatingIterableIterator consumer : iterators) {
				try {
					consumer.queue.put(n);
				} catch (InterruptedException e) {
					log.warn("Interrupted waiting on another iterator", e);
				}
			}
		}
	}
	private synchronized void removeIterator(DuplicatingIterableIterator consumer) {
		iterators.remove(consumer);
	}
	private class DuplicatingIterableIterator extends AbstractIterator<T> implements Closeable {
		private BlockingQueue<T> queue = new ArrayBlockingQueue<T>(bufferSize);
		private boolean isClosed = false;
		@Override
		protected T computeNext() {
			if (isClosed) {
				return endOfData();
			}
			T result = queue.poll();
			if (result != null) {
				return result;
			}
			consumeNext(this);
			result = queue.poll();
			if (result != null) {
				return result;
			}
			assert(queue.isEmpty());
			assert(!it.hasNext());
			return endOfData();
		}
		@Override
		public void close() {
			isClosed= true;
			removeIterator(this);
		}
	}
}
