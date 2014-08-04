package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicReference;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Wrapper iterator that uses a background thread to read from a given source iterator.
 * Unlike 
 * @author Daniel Cameron
 *
 */
public class AsyncBufferedIterator<T> implements CloseableIterator<T> {
	private static volatile int threadsCreated = 0;
	private final Log log = Log.getInstance(AsyncBufferedIterator.class);
    private final Thread reader;
    private final ReaderRunnable readerRunnable;
    private final AtomicReference<Throwable> ex = new AtomicReference<Throwable>(null);
    private final Iterator<T> underlying;
	private final BlockingQueue<List<Object>> buffer;
	private boolean closeCalled = false;
	private final int bufferSize;
    private PeekingIterator<Object> currentBuffer = Iterators.peekingIterator(ImmutableList.<Object>of().iterator());
	private static final Object eos = new Object(); // End of stream sentinel
	public AsyncBufferedIterator(Iterator<T> iterator, int bufferCount, int bufferSize) {
		if (iterator == null) throw new IllegalArgumentException();
		if (bufferCount <= 0 || bufferSize <= 0) throw new IllegalArgumentException("Non-zero buffer size required.");
		this.underlying = iterator;
		this.buffer = new ArrayBlockingQueue<List<Object>>(bufferCount);
		this.bufferSize = bufferSize;
        this.readerRunnable = new ReaderRunnable();
        this.reader = new Thread(readerRunnable, getThreadNamePrefix() + threadsCreated++);
        this.reader.setDaemon(true);
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
			reader.join();
		} catch (InterruptedException ie) { }
	}
	private void syncClose() {
		if (underlying instanceof Closeable) {
			try {
				((Closeable)underlying).close();
			} catch (Exception e) {
				log.error(e, "Error closing SAMRecord input stream");
			}
		}
	}
	@Override
	public boolean hasNext() {
		if (closeCalled) return false;
		throwOnCallingThread();
		if (!currentBuffer.hasNext()) {
			try {
				currentBuffer = Iterators.peekingIterator(buffer.take().iterator());
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
		}
		return currentBuffer.hasNext() && currentBuffer.peek() != eos;
	}
	@SuppressWarnings("unchecked")
	@Override
	public T next() {
		if (hasNext()) return (T)currentBuffer.next();
		throw new IllegalStateException("No more records");
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
        public void run() {
        	try {
        		while (underlying.hasNext()) {
		        	List<Object> readAhead = new ArrayList<Object>(bufferSize + 1);
		    		for (int i = 0; i < bufferSize; i++) {
		    			if (!underlying.hasNext()) break;
		    			readAhead.add(underlying.next());
		    		}
		    		if (!underlying.hasNext()) {
		    			readAhead.add(eos);
		    		}
		    		buffer.put(readAhead);
        		}
        	} catch (InterruptedException ie) {
        	} catch (Throwable t) {
        		ex.compareAndSet(null, t);
        	} finally {
        		syncClose(); 
        	}
        }
    }
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
