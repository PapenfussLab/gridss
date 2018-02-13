package au.edu.wehi.idsv.util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

/**
 * Wrapper iterator that uses a background thread to read from a given source iterator.
 * 
 * @author Daniel Cameron
 *
 */
public class AsyncBufferedIterator<T> implements CloseableIterator<T>, PeekingIterator<T>, AutoCloseable {
	private static AtomicInteger threadsCreated = new AtomicInteger(0);
	private static final Log log = Log.getInstance(AsyncBufferedIterator.class);
    private final Thread reader;
    private final ReaderRunnable readerRunnable;
    private final AtomicReference<Throwable> ex = new AtomicReference<Throwable>(null);
    private final Iterator<T> underlying;
	private final BlockingQueue<List<Object>> buffer;
	private boolean closeCalled = false;
	private final int batchSize;
    private PeekingIterator<Object> currentBuffer = Iterators.peekingIterator(ImmutableList.<Object>of().iterator());
	private static final Object eos = new Object(); // End of stream sentinel
	/**
	 * Creates a new iterator that traverses the given iterator on a background thread
	 * @param iterator iterator to traverse
	 * @param bufferCount number of read-ahead buffers
	 * @param batchSize size of each read-ahead buffer. A larger batch size will increase throughput and latency.
	 */
	public AsyncBufferedIterator(Iterator<T> iterator, int bufferCount, int batchSize) {
		this(iterator, null, bufferCount, batchSize);
	}
	public AsyncBufferedIterator(Iterator<T> iterator, String description) {
		this(iterator, description, gridss.Defaults.ASYNC_BUFFERS, gridss.Defaults.ASYNC_BUFFER_SIZE);
	}
	public AsyncBufferedIterator(Iterator<T> iterator, String description, int bufferCount, int batchSize) {
		if (iterator == null) throw new IllegalArgumentException();
		if (bufferCount <= 0 || batchSize <= 0) throw new IllegalArgumentException("Buffer size must be at least 1.");
		this.underlying = iterator;
		this.buffer = new ArrayBlockingQueue<List<Object>>(bufferCount);
		this.batchSize = batchSize;
        this.readerRunnable = new ReaderRunnable();
        this.reader = new Thread(readerRunnable, description == null ? getThreadNamePrefix() + threadsCreated.incrementAndGet() : description);
        this.reader.setDaemon(true);
        log.debug("Starting thread " + this.reader.getName());
        this.reader.start();
	}
	protected String getThreadNamePrefix() {
		return "AsyncBufferedIterator";
	}
	@Override
	public void close() {
		closeCalled = true;
		try {
			reader.interrupt();
			buffer.clear(); // flush buffer so EOS indicator can be written if writer is blocking
			reader.join();
		} catch (InterruptedException ie) { }
	}
	private void syncClose() {
		CloserUtil.close(underlying);
	}
	@Override
	public boolean hasNext() {
		throwOnCallingThread();
		if (closeCalled) return false;
		if (!currentBuffer.hasNext()) {
			try {
				currentBuffer = Iterators.peekingIterator(buffer.take().iterator());
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
			// rethrow any exceptions raised on the background thread while we were blocking on the next record
			throwOnCallingThread();
		}
		return currentBuffer.hasNext() && currentBuffer.peek() != eos;
	}
	@SuppressWarnings("unchecked")
	@Override
	public T next() {
		if (hasNext()) return (T)currentBuffer.next();
		throw new NoSuchElementException("next");
	}
	@SuppressWarnings("unchecked")
	@Override
	public T peek() {
		if (!hasNext()) throw new NoSuchElementException("peek");
		return (T)currentBuffer.peek();
	}
	private final void throwOnCallingThread() {
        final Throwable t = this.ex.get();
        if (t != null) {
            if (t instanceof Error) throw (Error) t;
            if (t instanceof RuntimeException) throw (RuntimeException) t;
            else throw new RuntimeException(t);
        }
    }
	/**
     * Reads the given iterator and passing back to the calling thread
     * in chunks
     */
    private class ReaderRunnable implements Runnable {
    	private boolean eosWritten = false;
        public void run() {
        	try {
        		while (underlying.hasNext()) {
		        	List<Object> readAhead = new ArrayList<Object>(batchSize + 1);
		    		for (int i = 0; i < batchSize; i++) {
		    			if (!underlying.hasNext()) break;
		    			readAhead.add(underlying.next());
		    		}
		    		if (!underlying.hasNext()) {
		    			readAhead.add(eos);
		    			eosWritten = true;
		    		}
		    		buffer.put(readAhead);
        		}
        	} catch (InterruptedException ie) {
        		// log.debug("Thread interrupt received - closing on background thread.");
        	} catch (Throwable t) {
        		// when using async I/O, htsjdk BlockCompressedInputStream
        		// RuntimeException-wrapped InterruptedExceptionon the calling thread
        		// (i.e., this one).
        		boolean causedByInterruptedException = false;
        		for (Throwable cause = t.getCause(); cause != null; cause = cause.getCause()) {
        			if (cause instanceof InterruptedException) {
        				causedByInterruptedException = true;
        				break;
        			}
        		}
        		if (!causedByInterruptedException) {
        			ex.set(t);
        			throw new RuntimeException(t);
        		}
        	} finally {
        		syncClose();
        		Thread.interrupted(); // clear thread interrupt flag so we can write the eos indicator if needed
        		try {
        			if (!eosWritten) {
        				buffer.put(ImmutableList.of(eos));
        			}
				} catch (InterruptedException e2) {
					log.warn("Thread interrupt received whilst writing end of stream indicator");
				}
        	}
        }
    }
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	protected String getBackgroundThreadName() {
		return this.reader.getName();
	}
}
