package au.edu.wehi.idsv.util;

import java.util.Iterator;

/**
 * Counts the number of elements that have been consumed
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class CountingIterator<T> implements Iterator<T> {
	private final Iterator<T> underlying;
	private int count = 0;
	public int emitted() {
		return count;
	}
	public CountingIterator(Iterator<T> it) {
		this.underlying = it;
	}
	@Override
	public boolean hasNext() {
		return underlying.hasNext();
	}

	@Override
	public T next() {
		T n = underlying.next();
		count++;
		return n;
	}

	@Override
	public void remove() {
		underlying.remove();
	}
}
