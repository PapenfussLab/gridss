package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.CloserUtil;

import java.io.Closeable;
import java.util.Iterator;

import com.google.common.collect.AbstractIterator;

/**
 * Iterator that automatically closes the underlying resources when the end of stream has been reached.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class AutoClosingIterator<T> extends AbstractIterator<T> implements Closeable {
	private final Iterator<T> underlying;
	private boolean closed = false;
	public AutoClosingIterator(Iterator<T> it) {
		this.underlying = it;
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
		}
	}
}
