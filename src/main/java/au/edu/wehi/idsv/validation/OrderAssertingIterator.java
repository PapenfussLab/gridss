package au.edu.wehi.idsv.validation;

import java.util.Comparator;
import java.util.Iterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

/**
 * Asserts that the iterator is in the specified order.
 * 
 * No checking is performed if assertions are not enabled. 
 * 
 * @author cameron.d
 *
 */
public class OrderAssertingIterator<T> implements CloseableIterator<T> {
	private static final Log log = Log.getInstance(OrderAssertingIterator.class);
	private final Iterator<? extends T> it;
	private final Comparator<? super T> comparator;
	private T last = null;
	private boolean closed = false;
	public OrderAssertingIterator(Iterator<? extends T> input, Comparator<? super T> comparator) {
		this.it = input;
		this.comparator = comparator;
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}
	@Override
	public T next() {
		T current = it.next();
		assert(checkSortOrder(current));
		return current; 
	}
	private boolean checkSortOrder(T current) {
		boolean success = true;
		if (last != null && current != null && comparator.compare(last, current) > 0) {
			success = false;
			log.error("Sanity check failure: iterator not sorted. " + last.toString() + " encountered before " + current.toString());
		}
		last = current;;
		return success;
	}
	@Override
	public void remove() {
		it.remove();
	}
	@Override
	public void close() {
		if (!closed) {
			CloserUtil.close(it);
			closed = true;
		}
	}
}
