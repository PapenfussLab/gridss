package au.edu.wehi.idsv.util;

import java.util.AbstractList;
import java.util.RandomAccess;

/**
 * Sliding Window List
 * 
 * Only the maxCapacity elements with the highest index are retained. All other elements are null
 * 
 * @author cameron.d
 *
 * @param <T>
 */
public class SlidingWindowList<E> extends AbstractList<E> implements RandomAccess {
	/**
	 * Circular array backing store
	 */
	private E[] buffer;
	/**
	 * Next free element of the buffer
	 */
	private int headIndex;
	public int getWindowSize() {
		return buffer.length;
	}
	@SuppressWarnings (value="unchecked")
	public SlidingWindowList(int windowSize)
	{
		if (windowSize <= 0) throw new IllegalArgumentException("Window size must be positive");
		buffer = (E[])(new Object[windowSize]);
		headIndex = -1;
	}
	@Override
	public E get(int index) {
		if (index <= headIndex - buffer.length) return expiredAccess(index);
		if (index > headIndex) throw new IndexOutOfBoundsException("Index: "+index+", Size: "+size());
		return buffer[index % buffer.length];
	}
	@Override
	public int size() {
		return headIndex + 1;
	}
	@Override
    public E set(int index, E e) {
		if (index <= headIndex - buffer.length) return expiredAccess(index);
		// clear out all values required for resizing the array
		// we can stop when we've got to the index of this new element, or we've cleared the entire buffer
		for (int i = headIndex + 1; i < index && i < headIndex + buffer.length + 1; i++) {
			E item = buffer[i % buffer.length];
			if (item != null) {
				buffer[i % buffer.length] = null;
				onWindowExited(item);
			}
		}
		buffer[index % buffer.length] = e;
		headIndex = index;
		onWindowEntered(e);
		return e;
    }
	/**
	 * Called whenever a non-null element is removed from the sliding window 
	 * @param item
	 */
	protected void onWindowExited(E item) {
	}
	/**
	 * Called whenever a non-null element is added to the sliding window 
	 * @param item
	 */
	protected void onWindowEntered(E item) {
	}
    @Override
    public void add(int index, E e) {
    	if (index <= headIndex - buffer.length) expiredAccess(index);
    	set(index, e);
    }
    private E expiredAccess(int index) {
    	// what to do when the caller does something to an expired index.
    	// either a) throw exception
    	// or b) null
    	// the latter is more compatible with the List<> interface
    	// throw new IndexOutOfBoundsException("Index: "+index+", Size: "+size() + ", Valid: " + startingIndex);
    	return null;
    }
}
