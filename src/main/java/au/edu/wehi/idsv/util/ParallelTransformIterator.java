package au.edu.wehi.idsv.util;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.PriorityQueue;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Executor;
import java.util.function.Function;

import com.google.common.collect.Ordering;

/**
 * Performs a given transformation operation over all elements of an iterator.
 * The transform is applied to multiple iterator elements in parallel with
 * the order of the resultant iteration unchanged.
 * 
 * This class is not thread-safe and access from multiple threads should
 * be synchronised.
 * 
 * @author Daniel Cameron
 *
 */
public class ParallelTransformIterator<T, U> implements Iterator<U> {
	private static class TransformResult<U> {
		public TransformResult(final long ordinal, final U result) {
			this.ordinal = ordinal;
			this.result = result;
		}
		public final long ordinal;
		public final U result;
		@SuppressWarnings("rawtypes")
		public static Ordering<TransformResult> byOrdinal = Ordering.natural().onResultOf((TransformResult tr) -> tr.ordinal);
	}
	private final Iterator<T> it;
	private final Function<T, U> f;
	private final int lookahead;
	private final Executor threadpool;
	private final ArrayBlockingQueue<TransformResult<U>> completed;
	private final PriorityQueue<TransformResult<U>> results;
	/**
	 * Number of records that have been read from the underlying iterator
	 * but not yet returned from this iterator
	 */
	private int dispatched = 0;
	/**
	 * Ordinal of last record returned from this iterator
	 */
	private long lastOrdinal = -1;
	
	/**
	 * Instantiates a new iterator
	 * @param it underlying iterator
	 * @param f transform function
	 * @param lookahead number of record to process in parallel
	 */
	public ParallelTransformIterator(final Iterator<T> it, final Function<T, U> f, final int lookahead, Executor threadpool) {
		this.it = it;
		this.f = f;
		this.lookahead = lookahead;
		this.completed = new ArrayBlockingQueue<TransformResult<U>>(lookahead);
		this.results = new PriorityQueue<ParallelTransformIterator.TransformResult<U>>(TransformResult.byOrdinal);
		this.threadpool = threadpool;
	}

	@Override
	public boolean hasNext() {
		return dispatched > 0 || it.hasNext();
	}

	@Override
	public U next() {
		if (!hasNext()) throw new NoSuchElementException();
		// dispatching here increases our latency as we're always going to have
		// lookahead record in our buffers, but it improves throughput as we're
		// not waiting until we have no records dispatched before requeuing. 
		dispatch();
		while (results.isEmpty() || results.peek().ordinal != lastOrdinal + 1) {
			TransformResult<U> record;
			try {
				record = completed.take();
			} catch (InterruptedException e) {
				// what's the correct interrupt handling mechanism?
				// if we swalloe then raise Thread.currentThread().interrupt();
				// after getting our record, we might end up blocking
				// since the worker thread doing the work might also have
				// been interrupted
				throw new RuntimeException(e);
			}
			results.add(record);
		}
		U result = results.poll().result;
		dispatched--;
		lastOrdinal++;
		dispatch();
		return result;
	}
	/**
	 * Dispatches records until we have lookahead records.
	 */
	private void dispatch() {
		while (dispatched < lookahead && it.hasNext()) {
			T record = it.next();
			dispatched++;
			dispatch(lastOrdinal + dispatched, record);
		}
	}
	private void dispatch(final long ordinal, final T record) {
		threadpool.execute(() -> {
			completed.add(new TransformResult<U>(ordinal, f.apply(record)));
		});
	}
}
