package au.edu.wehi.idsv.util;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMRecord;

import java.util.*;

/**
 * Takes a grouped iterator and ungroups records
 */
public class UngroupingIterator<T, U extends Iterable<T>> implements Iterator<T> {
    private final PeekingIterator<Iterable<T>> it;
    private Iterator<T> currentIt = null;

    public UngroupingIterator(Iterator<U> it) {
        this.it = Iterators.peekingIterator(it);
    }

    private void ensureNext() {
        while ((currentIt == null || !currentIt.hasNext()) && it.hasNext()) {
            Iterable<T> current = it.next();
            if (current != null) {
                currentIt = current.iterator();
            } else {
                currentIt = null;
            }
        }
    }

    @Override
    public boolean hasNext() {
        ensureNext();
        return currentIt != null && currentIt.hasNext();
    }

    @Override
    public T next() {
        ensureNext();
        if (currentIt == null || !currentIt.hasNext()) {
            throw new NoSuchElementException();
        }
        return currentIt.next();
    }
}
