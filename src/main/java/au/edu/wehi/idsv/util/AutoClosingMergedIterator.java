package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.io.Closeable;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

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
	private List<AutoClosingIterator<T>> stillOpen;
	private Iterator<? extends T> merged;
	private boolean closed = false;
	
	public AutoClosingMergedIterator(Iterable<? extends Iterator<? extends T>> iterators, Comparator<? super T> comparator) {
		stillOpen = Lists.newArrayList(Iterables.transform(iterators, new Function<Iterator<? extends T>, AutoClosingIterator<T>>() {
			@Override
			public AutoClosingIterator<T> apply(Iterator<? extends T> input) {
				return new AutoClosingIterator<>(input);
			}
		}));
		merged = Iterators.mergeSorted(stillOpen, comparator);
	}
	@Override
	public boolean hasNext() {
		return !closed && merged.hasNext();
		
	}
	@Override
	public T next() {
		T n = merged.next();
		tryclose(false);
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
				if (!it.hasNext() || forceClosed) {
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
