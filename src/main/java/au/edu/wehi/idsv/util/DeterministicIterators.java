package au.edu.wehi.idsv.util;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.UnmodifiableIterator;

import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DeterministicIterators {
    public static <T> Iterator<T> mergeSorted(Iterable<? extends Iterator<? extends T>> iterators, Comparator<? super T> comparator) {
        List<Iterator<? extends T>> lit = Lists.newArrayList(iterators);
        List<Iterator<Record<T>>> rits = IntStream.range(0, lit.size()).mapToObj(i -> new RecordPackingIterator<T>(i, lit.get(i))).collect(Collectors.toList());
        UnmodifiableIterator<Record<T>> mergedit = Iterators.mergeSorted(rits, recordComparator(comparator));
        return new RecordUnpackingIterator<T>(mergedit);
    }
    private static class RecordPackingIterator<T> implements Iterator<Record<T>> {
        private final int iteratorIndex;
        private final Iterator<? extends T> underlying;
        public RecordPackingIterator(int iteratorIndex, Iterator<? extends T> underlying) {
            this.iteratorIndex = iteratorIndex;
            this.underlying = underlying;
        }
        @Override
        public boolean hasNext() {
            return underlying.hasNext();
        }

        @Override
        public Record<T> next() {
            return new Record<T>(iteratorIndex, underlying.next());
        }
    }
    private static class RecordUnpackingIterator<T> implements Iterator<T> {
        private final Iterator<? extends Record<? extends T>> underlying;
        public RecordUnpackingIterator(Iterator<? extends Record<? extends T>> underlying) {
            this.underlying = underlying;
        }
        @Override
        public boolean hasNext() {
            return underlying.hasNext();
        }

        @Override
        public T next() {
            return underlying.next().record;
        }
    }
    private static class Record<T> {
        public final int iteratorIndex;
        public final T record;

        private Record(int iteratorIndex, T record) {
            this.iteratorIndex = iteratorIndex;
            this.record = record;
        }
    }
    private static <T> Comparator<Record<T>> recordComparator(Comparator<? super T> comparator) {
        return (o1, o2) -> {
            int cmp = comparator.compare(o1.record, o2.record);
            if (cmp != 0) return cmp;
            return Integer.compare(o1.iteratorIndex, o2.iteratorIndex);
        };
    }
}
