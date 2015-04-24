package au.edu.wehi.idsv.util;

import java.util.AbstractList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * Creates a collection of collections that usually contain a single element
 * @author cameron.d
 *
 * @param <T>
 */
public class FirstOverflowList<T> extends AbstractList<List<T>> {
	private List<? extends T> primary;
	private List<? extends List<? extends T>> overflow;
	/**
	 * Creates a new collection
	 * @param primary collection contains first element of each collection. Cannot be null
	 * @param overflow collection containing additional elements of each collection. Can be null, and child collections can be null.
	 */
	public FirstOverflowList(List<? extends T> primary, List<? extends List<? extends T>> overflow) {
		if (primary == null) throw new IllegalArgumentException("primary cannot be null");
		assert(primary != null);
		this.primary = primary;
		this.overflow = overflow;
	}
	@Override
	public List<T> get(int index) {
		if (index < 0 || index >= size()) throw new IndexOutOfBoundsException();
		return new ChildList(index);
	}
	public T get(int index, int offset) {
		if (index < 0 || index >= size()) throw new IndexOutOfBoundsException();
		if (offset < 0 || offset >= size(index)) throw new IndexOutOfBoundsException();
		if (offset == 0) return primary.get(index);
		return overflow.get(index).get(offset - 1);
	}
	public Iterable<T> asFlattenedIterable() {
		return new FlattenedIterable();
	}
	@Override
	public int size() {
		return primary.size();
	}
	public int size(int offset) {
		if (overflow == null) return 1;
		Collection<? extends T> c = overflow.get(offset);
		if (c == null) return 1;
		return 1 + c.size();
	}
	private class ChildList extends AbstractList<T> {
		private int offset;
		public ChildList(int offset) {
			this.offset = offset;
		}
		@Override
		public int size() {
			return FirstOverflowList.this.size(offset);
		}
		@Override
		public T get(int index) {
			return FirstOverflowList.this.get(offset, index);
		}
	}
	public int linearIndex(int index, int indexIndex) {
		int pos = 0;
		for (int i = 0; i < index; i++) {
			pos += size(i);
		}
		pos += indexIndex;
		return pos;
	}
	private class FlattenedIterator implements Iterator<T> {
		private int offset;
		private int offsetOffset;
		public FlattenedIterator(int offset, int offsetOffset) {
			this.offset = offset;
			this.offsetOffset = offsetOffset;
		}
		@Override
		public boolean hasNext() {
			return offset < primary.size() && offsetOffset < FirstOverflowList.this.size(offset);
		}
		@Override
		public T next() {
			T result = FirstOverflowList.this.get(offset, offsetOffset);
			moveForward();
			return result;
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException("NYI");
		}
		private void moveForward() {
			offsetOffset++;
			if (FirstOverflowList.this.size(offset) <= offsetOffset) {
				offset++;
				offsetOffset = 0;
			}
		}
	}
	private class FlattenedIterable implements Iterable<T> {
		@Override
		public Iterator<T> iterator() {
			return new FlattenedIterator(0, 0);
		}
	}
}
