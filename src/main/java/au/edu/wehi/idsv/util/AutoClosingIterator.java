package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.ImmutableList;

/**
 * Iterator that automatically closes the underlying resources when the end of stream has been reached.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class AutoClosingIterator<T> implements Closeable, CloseableIterator<T> {
	private static final Log log = Log.getInstance(AutoClosingIterator.class);
	private Iterator<? extends T> underlying;
	private List<? extends Closeable> alsoClose;
	private boolean closed = false;
	public AutoClosingIterator(Iterator<? extends T> it) {
		this(it, ImmutableList.<Closeable>of());
	}
	public AutoClosingIterator(Iterator<? extends T> it, Iterable<? extends Closeable> alsoClose) {
		this.underlying = it;
		this.alsoClose = ImmutableList.copyOf(alsoClose);
	}
	@Override
	public void close() {
		try {
			if (!closed) {
				closed = true;
				CloserUtil.close(underlying);
				for (Closeable c : alsoClose) {
					CloserUtil.close(c);
				}
				underlying = null;
				alsoClose = null;
			}
		} catch (Exception e) {
			log.warn(e, "Error closing resource.");
		}
	}
	@Override
	public boolean hasNext() {
		return !closed && underlying != null && underlying.hasNext();
	}
	@Override
	public T next() {
		if (closed) throw new IllegalStateException("underlying iterator closed");
		return underlying.next();
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
