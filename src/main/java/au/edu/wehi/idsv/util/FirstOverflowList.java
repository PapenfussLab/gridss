package au.edu.wehi.idsv.util;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.lang.IndexOutOfBoundsException;

/**
 * Creates a collection of collections that usually contain a single element
 * @author cameron.d
 *
 * @param <T>
 */
public class FirstOverflowList<T> implements List<List<T>> {
	private List<? extends T> primary;
	private List<? extends List<? extends T>> overflow;
	/**
	 * Creates a new collection
	 * @param primary collection contains first element of each collection. Cannot be null
	 * @param overflow collection containing additional elements of each collection. Can be null, and child collections can be null.
	 */
	public FirstOverflowList(List<? extends T> primary, List<? extends List<? extends T>> overflow) {
		assert(primary != null);
		this.primary = primary;
		this.overflow = overflow;
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
	@Override
	public boolean isEmpty() {
		return primary.isEmpty();
	}
	@Override
	public boolean contains(Object o) {
		if (primary.contains(o)) return true;
		if (overflow != null) {
			for (Collection<? extends T> c : overflow) {
				if (c.contains(o)) return true;
			}
		}
		return false;
	}
	@Override
	public Iterator<List<T>> iterator() {
		return new FirstOverflowListIterator();
	}
	public Iterator<T> flattenedIterator() {
		return new FirstOverflowIterator();
	}
	@Override
	public Object[] toArray() {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public <U> U[] toArray(U[] a) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public boolean remove(Object o) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public boolean containsAll(Collection<?> c) {
		for (Object o : c) {
			if (!contains(o)) return false;
		}
		return true;
	}
	@Override
	public boolean removeAll(Collection<?> c) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public boolean retainAll(Collection<?> c) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public void clear() {
		throw new UnsupportedOperationException("NYI");
	}
	private class FirstOverflowListIterator implements Iterator<List<T>> {
		private int i = 0;
		@Override
		public boolean hasNext() {
			return i < primary.size();
		}
		@Override
		public List<T> next() {
			return new OffsetList(i++);
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException("NYI");
		}
	}
	private class OffsetList implements List<T> {
		private int offset;
		public OffsetList(int offset) {
			this.offset = offset;
		}
		@Override
		public int size() {
			return FirstOverflowList.this.size(offset);
		}
		@Override
		public boolean isEmpty() {
			return false;
		}
		@Override
		public boolean contains(Object o) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public Iterator<T> iterator() {
			return new FirstOverflowIterator(offset);
		}
		@Override
		public Object[] toArray() {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public <U> U[] toArray(U[] a) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public boolean add(T e) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public boolean remove(Object o) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public boolean containsAll(Collection<?> c) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public boolean addAll(Collection<? extends T> c) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public boolean addAll(int index, Collection<? extends T> c) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public boolean removeAll(Collection<?> c) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public boolean retainAll(Collection<?> c) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public void clear() {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public T get(int index) {
			return FirstOverflowList.this.get(offset, index);
		}
		@Override
		public T set(int index, T element) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public void add(int index, T element) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public T remove(int index) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public int indexOf(Object o) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public int lastIndexOf(Object o) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public java.util.ListIterator<T> listIterator() {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public java.util.ListIterator<T> listIterator(int index) {
			throw new UnsupportedOperationException("NYI");
		}
		@Override
		public List<T> subList(int fromIndex, int toIndex) {
			throw new UnsupportedOperationException("NYI");
		}
	}
	private class FirstOverflowIterator implements Iterator<T> {
		private int offset = 0;
		private int offsetOffset = 0;
		private boolean traverse = true;
		public FirstOverflowIterator() {
		}
		public FirstOverflowIterator(int offset) {
			this.offset = offset;
			traverse = false;
		}
		@Override
		public boolean hasNext() {
			return offset < primary.size() && offsetOffset < FirstOverflowList.this.size(offset);
		}
		@Override
		public T next() {
			T result = FirstOverflowList.this.get(offset, offsetOffset);
			offsetOffset++;
			if (traverse && FirstOverflowList.this.size(offset) <= offsetOffset) {
				offset++;
				offsetOffset = 0;
			}
			return result;
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException("NYI");
		}
	}
	private class FirstOverflowFlattenedIterable implements Iterable<T> {
		@Override
		public Iterator<T> iterator() {
			return new FirstOverflowIterator();
		}
	}
	@Override
	public boolean add(List<T> e) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public boolean addAll(Collection<? extends List<T>> c) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public boolean addAll(int index, Collection<? extends List<T>> c) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public List<T> get(int index) {
		return new OffsetList(index);
	}
	public T get(int index, int offset) {
		if (index < 0 || index >= size()) throw new IndexOutOfBoundsException();
		if (offset < 0 || offset >= size(index)) throw new IndexOutOfBoundsException();
		if (offset == 0) return primary.get(index);
		return overflow.get(index).get(offset - 1);
	}
	@Override
	public List<T> set(int index, List<T> element) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public void add(int index, List<T> element) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public List<T> remove(int index) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public int indexOf(Object o) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public int lastIndexOf(Object o) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public ListIterator<List<T>> listIterator() {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public ListIterator<List<T>> listIterator(int index) {
		throw new UnsupportedOperationException("NYI");
	}
	@Override
	public List<List<T>> subList(int fromIndex, int toIndex) {
		throw new UnsupportedOperationException("NYI");
	}
	public Iterable<T> asFlattenedIterable() {
		return new FirstOverflowFlattenedIterable();
	}
}
