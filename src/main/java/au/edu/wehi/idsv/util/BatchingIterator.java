package au.edu.wehi.idsv.util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Groups records in the underlying iterator into fixed size batches
 * @param <T>
 */
public class BatchingIterator<T> implements Iterator<List<T>> {
    private final Iterator<T> underlying;
    private final int batchSize;

    public BatchingIterator(Iterator<T> underlying, int batchSize) {
        if (batchSize < 1 ) throw new IllegalArgumentException("batchSize must be positive");
        this.underlying = underlying;
        this.batchSize = batchSize;
    }

    @Override
    public boolean hasNext() {
        return underlying.hasNext();
    }

    @Override
    public List<T> next() {
        List<T> list = new ArrayList<>(batchSize);
        for (int i = 0; i < batchSize && underlying.hasNext(); i++) {
            list.add(underlying.next());
        }
        return list;
    }
}
