package au.edu.wehi.idsv.util;

import java.util.ArrayDeque;
import java.util.Iterator;

/**
 * Buffers the given iterator.
 * 
 * This class is a workaround for https://github.com/samtools/htsjdk/issues/760 and can be deleted
 * when this issue is resolved.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class BufferedIterator<T> implements Iterator<T> {
	private final Iterator<T> underlying;
	private final int bufferSize;
	private final ArrayDeque<T> buffer;
	public BufferedIterator(final Iterator<T> underlying, final int bufferSize) {
		this.underlying = underlying;
		this.bufferSize = bufferSize;
		this.buffer = new ArrayDeque<>(bufferSize);
		fillBuffer();
	}
	private void fillBuffer() {
		while (underlying.hasNext() && buffer.size() < bufferSize) {
			buffer.addLast(underlying.next());
		}
	}
	@Override
	public boolean hasNext() {
		return !buffer.isEmpty();
	}
	@Override
	public T next() {
		T r = buffer.removeFirst();
		fillBuffer();
		return r;
	}
}
