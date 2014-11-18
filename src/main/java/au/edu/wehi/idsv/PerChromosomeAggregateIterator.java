package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;
import java.util.NoSuchElementException;

import com.google.common.base.Function;

/**
 * Lazy iterator that iterates over all per chromosome evidence
 * only opening new file when required and closing as soon as possible 
 * @author Daniel Cameron
 *
 */
public class PerChromosomeAggregateIterator<T> implements CloseableIterator<T> {
	private final SAMSequenceDictionary dict;
	private final Function<String, Iterator<T>> perChrIterator;
	public PerChromosomeAggregateIterator(SAMSequenceDictionary dict, Function<String, Iterator<T>> function) {
		this.dict = dict;
		this.perChrIterator = function;
	}
	private int currentReferenceIndex = -1;
	private Iterator<T> currentSource = null;
	private boolean hasMoreContigs() {
		return currentReferenceIndex < dict.size() - 1;
	}
	private boolean currentHasNext() {
		if (currentSource != null) {
			if (currentSource.hasNext()) {
				return true;
			} else {
				// close as soon as we know there's no more records left
				closeCurrent();
			}
		}
		return false;
	}
	private boolean advanceContig() {
		assert(hasMoreContigs());
		closeCurrent();
		currentReferenceIndex++;
		currentSource = perChrIterator.apply(dict.getSequence(currentReferenceIndex).getSequenceName());
		return true;
	}
	private void closeCurrent() {
		if (currentSource != null) {
			CloserUtil.close(currentSource);
		}
		currentSource = null;
	}
	@Override
	public boolean hasNext() {
		while (!currentHasNext() && hasMoreContigs()) {
			advanceContig();
		}
		return currentHasNext();
	}
	@Override
	public T next() {
		if (hasNext()) return currentSource.next();
		throw new NoSuchElementException();
	}
	@Override
	public void close() {
		closeCurrent();
		currentReferenceIndex = dict.size();
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}