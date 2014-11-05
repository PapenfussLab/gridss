package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.io.Closeable;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

/**
 * Iterator that automatically closes the the underlying resources when their respective end of stream has been reached.
 * 
 * @author Daniel Cameron
 *
 * @param <T>
 */
public class AutoClosingMergedIterator<T> extends AbstractIterator<T> implements Closeable, CloseableIterator<T> {
	private List<Iterator<? extends T>> stillOpen;
	private Iterator<? extends T> merged;
	
	public AutoClosingMergedIterator(Iterable<? extends Iterator<? extends T>> iterators, Comparator<? super T> comparator) {
		stillOpen = Lists.newArrayList(iterators);
		merged = Iterators.mergeSorted(stillOpen, comparator);
	}
	@Override
	protected T computeNext() {
		if (merged == null || !merged.hasNext()) {
			tryclose(true);
			return endOfData();
		}
		tryclose(false);
		return merged.next();
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
				Iterator<? extends T> it = stillOpen.get(i);
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
}
