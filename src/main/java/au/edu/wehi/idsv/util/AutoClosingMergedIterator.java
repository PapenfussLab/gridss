package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.validation.OrderAssertingIterator;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

/**
 * Iterator that automatically closes the the underlying resources when their respective end of stream has been reached.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class AutoClosingMergedIterator<T> implements Closeable, CloseableIterator<T> {
	private Comparator<? super T> comparator;
	private List<AutoClosingIterator<T>> stillOpen;
	private Iterator<? extends T> merged;
	private boolean closed = false;
	private T lastEmitted = null;
	public AutoClosingMergedIterator(final Iterable<? extends Iterator<? extends T>> iterators, final Comparator<? super T> comparator) {
		this.comparator = comparator;
		this.stillOpen = Lists.newArrayList(Iterables.transform(iterators, new Function<Iterator<? extends T>, AutoClosingIterator<T>>() {
			@Override
			public AutoClosingIterator<T> apply(Iterator<? extends T> input) {
				List<Closeable> toClose = new ArrayList<Closeable>();
				if (input instanceof Closeable) {
					toClose.add((Closeable)input);
				}
				return new AutoClosingIterator<T>(new OrderAssertingIterator<T>(input, comparator), toClose);
			}
		}));
		this.merged = Iterators.mergeSorted(stillOpen, comparator);
	}
	@Override
	public boolean hasNext() {
		return !closed && merged.hasNext();
		
	}
	@Override
	public T next() {
		T n = merged.next();
		tryclose(false);
		if (lastEmitted != null && comparator.compare(lastEmitted, n) > 0) {
			throw new IllegalStateException(String.format("Unable to merge out of order sequences. %s emitted before %s", lastEmitted, n));
		}
		lastEmitted = n;
		return n;
	}
	@Override
	public void close() {
		tryclose(true);
	}
	public void tryclose(boolean forceClosed) {
		if (merged == null) {
			return;
		} else if (forceClosed) {
			merged = null;
		}
		if (stillOpen != null) {
			for (int i = stillOpen.size() - 1; i >= 0 ; i--) {
				AutoClosingIterator<T> it = stillOpen.get(i);
				if (forceClosed || !it.hasNext()) {
					CloserUtil.close(it);
					stillOpen.remove(i);
				}
			}
			if (stillOpen.size() == 0) {
				stillOpen = null;
			}
		}
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
