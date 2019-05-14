package au.edu.wehi.idsv.util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.PeekingIterator;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;

/**
 * Duplicates the given iterator, feeding internal buffers from a background thread
 * 
 * This wrapper is thread-safe.
 * 
 * <b>Separate consumer threads are required as
 * iterator calls block the calling thread when sufficiently
 * far ahead of other iterators.
 * </b>
 * @author Daniel Cameron
 *
 */
public class DuplicatingIterable<T> implements Iterable<T> {
	private static final Log log = Log.getInstance(DuplicatingIterable.class);
	private static final Object endofstream = new Object();
	private static final AtomicInteger threadCount = new AtomicInteger(0);
	private final Iterator<T> it;
	private final List<DuplicatingIterableIterator> iterators = new ArrayList<DuplicatingIterableIterator>();
	private final List<BlockingQueue<Object>> queues = new ArrayList<BlockingQueue<Object>>();
	private int iteratorsRequested = 0;
	private FeedingThread thread;
	private volatile Exception error = null;
	
	/**
	 * Duplicates an iterator
	 * @param nIterators number of consuming iterators
	 * @param it underlying iterator
	 * @param maxIteratorDifference maximum number of records an iterator can traverse before blocking to wait
	 * for other iterators to catch up
	 */
	public DuplicatingIterable(int nIterators, Iterator<T> it,  int maxIteratorDifference) {
		if (it == null) throw new IllegalArgumentException();
		if (maxIteratorDifference <= 0) throw new IllegalArgumentException("buffer size must be greater than zero.");
		this.it = it;
		for (int i = 0; i < nIterators; i++) {
			queues.add(new ArrayBlockingQueue<Object>(maxIteratorDifference));
			iterators.add(new DuplicatingIterableIterator(queues.get(i)));
		}
		this.thread = new FeedingThread();
		this.thread.setName(String.format("DuplicatingIterable-%d", threadCount.incrementAndGet()));
		this.thread.start();
	}
	/**
	 * Creates a new iterator
	 */
	@Override
	public synchronized PeekingIterator<T> iterator() {
		if (iteratorsRequested >= iterators.size()) throw new IllegalStateException(String.format("Already created %d iterators", iterators.size()));
		return iterators.get(iteratorsRequested++);
	}
	private class FeedingThread extends Thread {
		@Override
		public void run() {
			try {
				while (it.hasNext()) {
					T n = it.next();
					for (BlockingQueue<Object> queue : queues) {
						queue.put(n);
					}
				}
				eos();
			} catch (InterruptedException e) {
				log.warn("Interrupted waiting to feed next record - ending stream early");
				for (BlockingQueue<Object> queue : queues) {
					queue.clear();
					try {
						eos();
					} catch (InterruptedException e1) {
						log.error("Sanity check failure: end of stream writing should not have blocked.");
					}
				}
			} catch (Exception e) {
				log.error("Error traversing iterator", e);
				error = e;
				try {
					eos();
				} catch (InterruptedException e1) {
				}
			}
		}
		private void eos() throws InterruptedException {
			for (BlockingQueue<Object> queue : queues) {
				queue.put(endofstream);
			}
		}
	}
	private class DuplicatingIterableIterator implements PeekingIterator<T> {
		private final BlockingQueue<Object> queue;
		/**
		 * Since BlockingQueue does not allow nulls, we can use it as a
		 * sentinal as to whether we have cached the next result 
		 */
		private Object nextRecord = null;
		public DuplicatingIterableIterator(BlockingQueue<Object> queue) {
			this.queue = queue;
		}
		private void ensureNext() {
			if (nextRecord == endofstream) return;
			if (nextRecord == null) {
				try {
					nextRecord = queue.take();
				} catch (InterruptedException e) {
					log.debug("Interrupted waiting for next record");
					throw new RuntimeException(e);
				}
			}
			if (error != null) {
				throw new RuntimeException(error);
			}
		}
		@Override
		public boolean hasNext() {
			ensureNext();
			return nextRecord != endofstream;
		}
		@SuppressWarnings("unchecked")
		@Override
		public T peek() {
			if (!hasNext()) throw new NoSuchElementException();
			return (T)nextRecord;
		}
		@SuppressWarnings("unchecked")
		@Override
		public T next() {
			if (!hasNext()) throw new NoSuchElementException();
			ensureNext();
			T result = (T)nextRecord;
			nextRecord = null; // invalidate cached record
			return result;
		}
		@Override
		public void remove() {
			throw new IllegalStateException();
		}
	}
}
