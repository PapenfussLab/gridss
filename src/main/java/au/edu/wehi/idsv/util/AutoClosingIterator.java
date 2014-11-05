package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.io.Closeable;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableList;

/**
 * Iterator that automatically closes the underlying resources when the end of stream has been reached.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class AutoClosingIterator<T> extends AbstractIterator<T> implements Closeable, CloseableIterator<T> {
	private final Iterator<? extends T> underlying;
	private final List<Closeable> alsoClose;
	private boolean closed = false;
	
	public AutoClosingIterator(Iterator<? extends T> it) {
		this.underlying = it;
		this.alsoClose = ImmutableList.of();
		if (!underlying.hasNext()) close();
	}
	public AutoClosingIterator(Iterator<T> it, Iterable<Closeable> alsoClose) {
		this.underlying = it;
		this.alsoClose = ImmutableList.copyOf(alsoClose);
		if (!underlying.hasNext()) close();
	}
	@Override
	protected T computeNext() {
		if (closed || !underlying.hasNext()) {
			close();
			return endOfData();
		}
		return underlying.next();
	}
	@Override
	public void close() {
		if (!closed) {
			closed = true;
			CloserUtil.close(underlying);
			for (Closeable c : alsoClose) {
				CloserUtil.close(c);
			}
		}
	}
}
