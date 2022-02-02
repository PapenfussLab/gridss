package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.util.IntervalUtil;
import com.google.common.collect.ImmutableList;
import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectSortedMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

/**
 * Kmer interval lookup in which no records overlap
 */
public abstract class KmerIntervalLookup<T> {
    /**
     * Lookup of node starts. Secondary key is end()
     * Values are the KmerNode itself (when there is only 1) or a Int2ObjectSortedMap
     */
    private final Long2ObjectOpenHashMap<Object> kmerLookup = new Long2ObjectOpenHashMap<>();
    private int size = 0;
    protected abstract int getStart(T node);
    protected abstract int getEnd(T node);
    protected abstract long getKmer(T node);
    public void add(T node) {
        long kmer = getKmer(node);
        Object x = kmerLookup.get(kmer);
        Int2ObjectSortedMap<T> positionLookup;
        if (x == null) {
            kmerLookup.put(kmer, node);
            return;
        } else if (x instanceof Int2ObjectSortedMap) {
            positionLookup = (Int2ObjectSortedMap<T>) x;
        } else {
            T existing = (T)x;
            positionLookup = new Int2ObjectRBTreeMap<>();
            kmerLookup.put(kmer, positionLookup);
            positionLookup.put(getEnd(existing), existing);
        }
        positionLookup.put(getEnd(node), node);
    }
    public void remove(T node) {
        long kmer = getKmer(node);
        Object x = kmerLookup.get(kmer);
        if (x instanceof Int2ObjectSortedMap) {
            Int2ObjectSortedMap<T> positionLookup = (Int2ObjectSortedMap<T>) x;
            T found = positionLookup.remove(getEnd(node));
            assert (found != null);
            if (positionLookup.isEmpty()) {
                kmerLookup.remove(kmer);
            }
        } else {
            Object found = kmerLookup.remove(kmer);
            assert(found != null);
        }
    }
    /**
     * Gets the KmerNode that overlaps exactly
     * @param kmer
     * @param start
     * @param end
     * @return
     */
    public T get(long kmer, int start, int end) {
        Object x = kmerLookup.get(kmer);
        T node = null;
        if (x instanceof Int2ObjectSortedMap) {
            Int2ObjectSortedMap<T> positionLookup = (Int2ObjectSortedMap<T>) x;
            node = positionLookup.get(end);
        } else {
            node = (T)x;
        }
        if (node != null && (getStart(node) != start || getEnd(node) != end)) {
            // doesn't overlap exactly
            return null;
        }
        return node;
    }
    public List<T> getOverlapping(long kmer, int start, int end) {
        Object x = kmerLookup.get(kmer);
        T node = null;
        if (x == null) {
            return Collections.EMPTY_LIST;
        } if (x instanceof Int2ObjectSortedMap) {
            Int2ObjectSortedMap<T> positionLookup = (Int2ObjectSortedMap<T>) x;
            positionLookup = positionLookup.tailMap(start);
            ArrayList result = new ArrayList();
            Iterator<T> it = positionLookup.values().stream().iterator();
            while (it.hasNext()) {
                node = it.next();
                if (IntervalUtil.overlapsClosed(start, end, getStart(node), getEnd(node))) {
                    result.add(node);
                } else {
                    break;
                }
            }
            return result;
        } else {
            node = (T)x;
            if (IntervalUtil.overlapsClosed(start, end, getStart(node), getEnd(node))) {
                return ImmutableList.of(node);
            }
            return Collections.EMPTY_LIST;
        }
    }
    public Stream<T> stream() {
        return kmerLookup.values().stream().flatMap(o -> stream(o));
    }
    private Stream<T> stream(Object x) {
        if (x == null) return Stream.empty();
        if (x instanceof Int2ObjectSortedMap) return ((Int2ObjectSortedMap<T>)x).values().stream();
        return Stream.of((T)x);
    }
}
