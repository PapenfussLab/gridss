package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerNode;
import au.edu.wehi.idsv.util.IntervalUtil;
import com.google.common.collect.Iterables;
import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectSortedMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * KmerNode lookup in which no records overlap
 */
public class KmerNodeNonOverlappingLookup<T extends KmerNode> {
    protected final int k;
    /**
     * Lookup of node starts. Secondary key is firstEnd()
     * Values are the KmerNode itself (when there is only 1) or a Int2ObjectSortedMap
     */
    private final Long2ObjectOpenHashMap<Object> startKmerLookup = new Long2ObjectOpenHashMap<>();
    /**
     * Lookup of node ends. Secondary key is -lastStart()
     * Values are the KmerNode itself (when there is only 1) or a Int2ObjectSortedMap
     * lastStart() is negated so our traversal order is backwards when we iterate over values()
     */
    private final Long2ObjectOpenHashMap<Object> endKmerLookup = new Long2ObjectOpenHashMap<>();
    private int size = 0;
    public KmerNodeNonOverlappingLookup(int k) {
        this.k = k;
    }

    public List<T> prevNodes(T right) {
        return prevNodes(right, false);
    }
    public List<T> prevNodes(T right, boolean abortIfNotSingleUniqueFullWidth) {
        List<T> adj = new ArrayList<T>(4);
        for (long kmer : KmerEncodingHelper.prevStates(k, right.firstKmer())) {
            Object lookup = endKmerLookup.get(kmer);
            if (lookup == null) {
                // nothing to do
            } else if (lookup instanceof Int2ObjectSortedMap) {
                Int2ObjectSortedMap<T> candidates = ((Int2ObjectSortedMap<T>)lookup).tailMap(-(right.firstEnd() - 1));
                for (T left : candidates.values()) {
                    if (left.lastEnd() < right.firstStart() - 1) {
                        break;
                    } else {
                        assert(IntervalUtil.overlapsClosed(left.lastStart() + 1, left.lastEnd() + 1, right.firstStart(), right.firstEnd()));
                        assert(KmerEncodingHelper.isNext(k, left.lastKmer(), right.firstKmer()));
                        if (abortIfNotSingleUniqueFullWidth && (adj.size() != 0 || !(left.lastStart() + 1 == right.firstStart() && left.lastEnd() + 1 == right.firstEnd()))) {
                            return null;
                        }
                        adj.add(left);
                    }
                }
            } else {
                T left = (T)lookup;
                if (IntervalUtil.overlapsClosed(left.lastStart() + 1, left.lastEnd() + 1, right.firstStart(), right.firstEnd())) {
                    if (abortIfNotSingleUniqueFullWidth && (adj.size() != 0 || !(left.lastStart() + 1 == right.firstStart() && left.lastEnd() + 1 == right.firstEnd()))) {
                        return null;
                    }
                    adj.add(left);
                }
            }
        }
        return adj;
    }
    public List<T> nextNodes(T left) {
        return nextNodes(left, false);
    }

    private List<T> nextNodes(T left, boolean abortIfNotSingleUniqueFullWidth) {
        List<T> adj = new ArrayList<T>(4);
        for (long kmer : KmerEncodingHelper.nextStates(k, left.lastKmer())) {
            Object lookup = startKmerLookup.get(kmer);
            if (lookup == null) {
                // nothing to do
            } else if (lookup instanceof Int2ObjectSortedMap) {
                Int2ObjectSortedMap<T> candidates = ((Int2ObjectSortedMap<T>)lookup).tailMap(left.lastStart() + 1);
                for (T right : candidates.values()) {
                    if (right.firstStart() > left.lastEnd() + 1) {
                        break;
                    } else {
                        assert(IntervalUtil.overlapsClosed(left.lastStart() + 1, left.lastEnd() + 1, right.firstStart(), right.firstEnd()));
                        assert(KmerEncodingHelper.isNext(k, left.lastKmer(), right.firstKmer()));
                        if (abortIfNotSingleUniqueFullWidth && (adj.size() != 0 || !(left.lastStart() + 1 == right.firstStart() && left.lastEnd() + 1 == right.firstEnd()))) {
                            return null;
                        }
                        adj.add(right);
                    }
                }
            } else {
                T right = (T)lookup;
                if (IntervalUtil.overlapsClosed(left.lastStart() + 1, left.lastEnd() + 1, right.firstStart(), right.firstEnd())) {
                    if (abortIfNotSingleUniqueFullWidth && (adj.size() != 0 || !(left.lastStart() + 1 == right.firstStart() && left.lastEnd() + 1 == right.firstEnd()))) {
                        return null;
                    }
                    adj.add(right);
                }
            }
        }
        return adj;
    }

    public T getUniqueFullWidthSuccessor(T left) {
        List<T> adj = nextNodes(left, true);
        return (adj != null && !adj.isEmpty()) ? adj.get(0) : null;
    }

    public T getUniqueFullWidthPredecessor(T right) {
        List<T> adj = prevNodes(right, true);
        return (adj != null && !adj.isEmpty()) ? adj.get(0) : null;
    }
    public void remove(T kn) {
        //sanityCheckContains(kn);
        removeFromStartKmerLookup(kn);
        removeFromEndKmerLookup(kn);
        size--;
        //sanityCheck();
    }
    public void removeFromStartKmerLookup(T kn) {
        long kmer = kn.firstKmer();
        Object x = startKmerLookup.get(kmer);
        if (x instanceof Int2ObjectSortedMap) {
            Int2ObjectSortedMap<T> positionLookup = (Int2ObjectSortedMap<T>) x;
            T found = positionLookup.remove(kn.firstEnd());
            assert (found != null);
            if (positionLookup.isEmpty()) {
                startKmerLookup.remove(kmer);
            }
        } else {
            Object found = startKmerLookup.remove(kmer);
            assert (found != null);
        }
    }
    public void removeFromEndKmerLookup(T kn) {
        long kmer = kn.lastKmer();
        Object x = endKmerLookup.get(kmer);
        if (x instanceof Int2ObjectSortedMap) {
            Int2ObjectSortedMap<T> positionLookup = (Int2ObjectSortedMap<T>) x;
            T found = positionLookup.remove(-kn.lastStart());
            assert (found != null);
            if (positionLookup.isEmpty()) {
                endKmerLookup.remove(kmer);
            }
        } else {
            Object found = endKmerLookup.remove(kmer);
            assert (found != null);
        }
    }
    public void add(T kn) {
        addOrReplaceStartKmerLookup(kn, kn);
        addOrReplaceEndKmerLookup(kn, kn);
        size++;
        //sanityCheck();
    }

    private void addOrReplaceStartKmerLookup(T key, T value) {
        long kmer = key.firstKmer();
        Object x = startKmerLookup.get(kmer);
        Int2ObjectSortedMap<T> positionLookup;
        if (x == null) {
            startKmerLookup.put(kmer, value);
            return;
        } else if (x instanceof Int2ObjectSortedMap) {
            positionLookup = (Int2ObjectSortedMap<T>) x;
        } else {
            T existing = (T)x;
            if (existing.firstEnd() == key.firstEnd()) {
                assert(existing != value);
                // We are doing a replace
                // since we only have 1 record, we can just overwrite it
                startKmerLookup.put(kmer, value);
                return;
            }
            positionLookup = new Int2ObjectRBTreeMap<>();
            startKmerLookup.put(kmer, positionLookup);
            positionLookup.put(existing.firstEnd(), existing);
        }
        positionLookup.put(key.firstEnd(), value);
    }

    private void addOrReplaceEndKmerLookup(T key, T value) {
        long kmer = key.lastKmer();
        Object x = endKmerLookup.get(kmer);
        Int2ObjectSortedMap<T> positionLookup;
        if (x == null) {
            endKmerLookup.put(kmer, value);
            return;
        } else if (x instanceof Int2ObjectSortedMap) {
            positionLookup = (Int2ObjectSortedMap<T>) x;
        } else {
            T existing = (T)x;
            if (existing.lastStart() == key.lastStart()) {
                assert(existing != value);
                // We are doing a replace
                // since we only have 1 record, we can just overwrite it
                endKmerLookup.put(kmer, value);
                return;
            }
            positionLookup = new Int2ObjectRBTreeMap<>();
            endKmerLookup.put(kmer, positionLookup);
            positionLookup.put(-existing.lastStart(), existing);
        }
        positionLookup.put(-key.lastStart(), value);
    }

    public void replace(T toRemove, T toAdd) {
        assert(toRemove.firstKmer() == toAdd.firstKmer());
        assert(toRemove.lastKmer() == toAdd.lastKmer());
        assert(toRemove.firstStart() == toAdd.firstStart());
        assert(toRemove.firstEnd() == toAdd.firstEnd());
        assert(toRemove.lastStart() == toAdd.lastStart());
        assert(toRemove.lastEnd() == toAdd.lastEnd());
        addOrReplaceStartKmerLookup(toRemove, toAdd);
        addOrReplaceEndKmerLookup(toRemove, toAdd);
    }

    public boolean isEmpty() {
        return startKmerLookup.isEmpty();
    }

    public int size() {
        return size;
    }

    /**
     * Updates the lookup based on the right node getting merged into the left
     *
     * WARNING: the lookup will be corrupt until left is actually updated to include right.
     * @param left node to be extended
     * @param right node to be replaced
     */
    public void adjustForMerge(T left, T right) {
        removeFromEndKmerLookup(left);
        removeFromStartKmerLookup(right);
        addOrReplaceEndKmerLookup(right, left); // this peforms a replace
        size--;
    }
    private void sanityCheckContains(KmerNode kn) {
        List<KmerNode> start = startKmerLookup.values()
                .stream()
                .flatMap(x -> (x instanceof Int2ObjectSortedMap) ? ((Int2ObjectSortedMap<T>)x).values().stream() : Stream.of((T)x))
                .filter(x -> x == kn)
                .collect(Collectors.toList());
        List<KmerNode> end = endKmerLookup.values()
                .stream()
                .flatMap(x -> (x instanceof Int2ObjectSortedMap) ? ((Int2ObjectSortedMap<T>)x).values().stream() : Stream.of((T)x))
                .filter(x -> x == kn)
                .collect(Collectors.toList());
        assert(start.size() == 1);
        assert(end.size() == 1);
    }
    public void sanityCheck() {
        List<KmerNode> start = startKmerLookup.values()
                .stream()
                .flatMap(x -> (x instanceof Int2ObjectSortedMap) ? ((Int2ObjectSortedMap<T>)x).values().stream() : Stream.of((T)x))
                .collect(Collectors.toList());
        List<KmerNode> end = endKmerLookup.values()
                .stream()
                .flatMap(x -> (x instanceof Int2ObjectSortedMap) ? ((Int2ObjectSortedMap<T>)x).values().stream() : Stream.of((T)x))
                .collect(Collectors.toList());
        assert(start.size() == end.size());
        assert(size == start.size());
        for (KmerNode kn : Iterables.concat(start, end)) {
            Object x = startKmerLookup.get(kn.firstKmer());
            assert(x != null);
            if (x instanceof Int2ObjectSortedMap) {
                assert(((Int2ObjectSortedMap<T>)x).get(kn.firstEnd()) == kn);
            } else {
                assert(x == kn);
            }
            x = endKmerLookup.get(kn.lastKmer());
            if (x instanceof Int2ObjectSortedMap) {
                assert(((Int2ObjectSortedMap<T>)x).get(-kn.lastStart()) == kn);
            } else {
                assert(x == kn);
            }
        }
    }
}
