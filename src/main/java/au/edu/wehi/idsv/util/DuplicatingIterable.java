package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.AbstractIterator;

/**
 * Duplicates the given iterator, feeding internal buffers from a background thread
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
	private static final Object endofstream = new Object();
	private static final AtomicInteger threadCount = new AtomicInteger(0);
	private final Iterator<T> it;
	private final List<DuplicatingIterableIterator> iterators = new ArrayList<DuplicatingIterableIterator>();
	private final List<BlockingQueue<Object>> queues = new ArrayList<BlockingQueue<Object>>();
	private int iteratorsRequested = 0;
	private FeedingThread thread;
	
	/**
	 * Duplicates an iterator
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
	 * Creates a new iterator from the current position
	 */
	@Override
	public synchronized Iterator<T> iterator() {
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
			}
		}
		private void eos() throws InterruptedException {
			for (BlockingQueue<Object> queue : queues) {
				queue.put(endofstream);
			}
		}
	}
	private class DuplicatingIterableIterator extends AbstractIterator<T> {
		private final BlockingQueue<Object> queue; 
		public DuplicatingIterableIterator(BlockingQueue<Object> queue) {
			this.queue = queue;
		}
		@SuppressWarnings("unchecked")
		@Override
		protected T computeNext() {
			Object o;
			try {
				o = queue.take();
				if (o != endofstream) return (T)o;
			} catch (InterruptedException e) {
				log.debug("Interrupted waiting for next record");
				throw new RuntimeException(e);
			}
			return endOfData();
		}
	}
}
