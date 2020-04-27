package au.edu.wehi.idsv.util;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * Takes a sorted iterator and groups records together
 */
public class GroupingIterator<T> implements Iterator<List<T>> {
    private final PeekingIterator<T> it;
    private final Comparator<T> comparator;

    public GroupingIterator(Iterator<T> it, Comparator<T> comparator) {
        this.it = Iterators.peekingIterator(it);
        this.comparator = comparator;
    }

    @Override
    public boolean hasNext() {
        return it.hasNext();
    }

    @Override
    public List<T> next() {
        T ref = it.next();
        List<T> records = new ArrayList<>();
        records.add(ref);
        while (it.hasNext() && comparator.compare(ref, it.peek()) == 0) {
            records.add(it.next());
        }
        return records;
    }
}
