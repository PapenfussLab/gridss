package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.NoSuchElementException;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * Filters the iterable to only results on the given chromsome
 * closing the underlying resources as soon as possible
 * 
 * @author Daniel Cameron
 *
 */
public class ChromosomeFilteringIterator<T extends DirectedEvidence> implements CloseableIterator<T> {
	private final Iterator<T> source;
	private final PeekingIterator<T> it;
	private final int referenceIndex;
	private boolean isSorted;
	private boolean closed = false;
	/**
	 * Creates an iterator that returns 
	 * @param it source iterator
	 * @param referenceIndex reference index to return
	 * @param isSorted input iterator is sorted by reference index
	 */
	public ChromosomeFilteringIterator(Iterator<T> it, int referenceIndex, boolean isSorted) {
		assert(referenceIndex >= 0);
		this.source = it;
		this.it = Iterators.peekingIterator(it);
		this.referenceIndex = referenceIndex;
		this.isSorted = isSorted;
	}
	@Override
	public boolean hasNext() {
		if (closed) return false;
		if (isSorted) {
			while (it.hasNext() && it.peek().getBreakendSummary().referenceIndex < referenceIndex) {
				// fast-forward advance to our chr
				it.next();
			}
			if (it.hasNext() && it.peek().getBreakendSummary().referenceIndex > referenceIndex) {
				// if we've nothing left on our chr, close immediately
				close();
				return false;
			}
		}
		while (it.hasNext() && it.peek().getBreakendSummary().referenceIndex != referenceIndex) {
			it.next();
		}
		if (!it.hasNext()) {
			close();
			return false;
		}
		return !closed && it.hasNext();
	}
	@Override
	public T next() {
		if (hasNext()) return it.next();
		throw new NoSuchElementException();
	}
	@Override
	public void close() {
		if (!closed) {
			CloserUtil.close(source);
		}
		closed = true;
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}