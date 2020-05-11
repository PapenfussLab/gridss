package au.edu.wehi.idsv.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Iterator equivalent of stream.flatMap()
 * @param <T> Type of iterator
 */
public class FlatMapIterator<T> implements Iterator<T> {
    private final Iterator<? extends Iterable<T>> underlying;
    private Iterator<T> currentIterator = Collections.emptyIterator();

    public FlatMapIterator(Iterator<? extends Iterable<T>> underlying) {
        this.underlying = underlying;
    }

    public void ensureCurrentIterator() {
        while (!currentIterator.hasNext() && underlying.hasNext()) {
            currentIterator = underlying.next().iterator();
        }
    }

    @Override
    public boolean hasNext() {
        ensureCurrentIterator();
        return currentIterator.hasNext();
    }

    @Override
    public T next() {
        ensureCurrentIterator();
        return currentIterator.next();
    }
}
