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
     * Lookup of node starts. Secondary key is the end of the start kmer position interval or the single value itself
     */
    private final Long2ObjectOpenHashMap<Object> startKmerLookup = new Long2ObjectOpenHashMap<>();
    /**
     * Lookup of node ends. Secondary key is the end of the end kmer position interval or the single value itself
     */
    private final Long2ObjectOpenHashMap<Object> endKmerLookup = new Long2ObjectOpenHashMap<>();
    private int size = 0;
    public KmerNodeNonOverlappingLookup(int k) {
        this.k = k;
    }

    public List<T> prevNodes(T node) {
        List<T> adj = new ArrayList<T>(4);
        for (long kmer : KmerEncodingHelper.prevStates(k, node.lastKmer())) {
            Object lookup = endKmerLookup.get(kmer);
            if (lookup == null) {
                // nothing to do
            } else if (lookup instanceof Int2ObjectSortedMap) {
                Int2ObjectSortedMap<T> candidates = ((Int2ObjectSortedMap<T>)lookup).tailMap(node.firstStart() - 1);
                for (T n : candidates.values()) {
                    if (n.firstStart() > node.firstEnd() - 1) {
                        break;
                    } else {
                        assert (IntervalUtil.overlapsClosed(n.lastStart() + 1, n.lastEnd() + 1, node.firstStart(), node.firstEnd()));
                        assert (KmerEncodingHelper.isNext(k, n.firstKmer(), node.lastKmer()));
                        adj.add(n);
                    }
                }
            } else {
                T n = (T)lookup;
                if (IntervalUtil.overlapsClosed(n.lastStart() + 1, n.lastEnd() + 1, node.firstStart(), node.firstEnd())) {
                    adj.add(n);
                }
            }
        }
        return adj;
    }
    public List<T> nextNodes(T node) {
        List<T> adj = new ArrayList<T>(4);
        for (long kmer : KmerEncodingHelper.nextStates(k, node.lastKmer())) {
            Object lookup = startKmerLookup.get(kmer);
            if (lookup == null) {
                // nothing to do
            } else if (lookup instanceof Int2ObjectSortedMap) {
                Int2ObjectSortedMap<T> candidates = ((Int2ObjectSortedMap<T>)lookup).tailMap(node.lastStart() + 1);
                for (T n : candidates.values()) {
                    if (n.lastStart() > node.lastEnd() + 1) {
                        break;
                    } else {
                        assert (IntervalUtil.overlapsClosed(node.lastStart() + 1, node.lastEnd() + 1, n.firstStart(), n.firstEnd()));
                        assert (KmerEncodingHelper.isNext(k, node.lastKmer(), n.firstKmer()));
                        adj.add(n);
                    }
                }
            } else {
                T n = (T)lookup;
                if (IntervalUtil.overlapsClosed(node.lastStart() + 1, node.lastEnd() + 1, n.firstStart(), n.firstEnd())) {
                    adj.add(n);
                }
            }
        }
        return adj;
    }

    public T getUniqueSuccessor(T node) {
        T candidate = null;
        for (long kmer : KmerEncodingHelper.nextStates(k, node.lastKmer())) {
            Object lookup = startKmerLookup.get(kmer);
            if (lookup == null) {
                // nothing to do
            } else if (lookup instanceof Int2ObjectSortedMap) {
                Int2ObjectSortedMap<T> candidates = ((Int2ObjectSortedMap<T>)lookup);
                T next = candidates.get(node.firstEnd() + 1);
                if (next != null && next.firstStart() + 1 == node.firstStart()) {
                    if (candidate != null) {
                        // found a second candidate
                        return null;
                    }
                    candidate = next;
                }
            } else {
                T next = (T)lookup;
                if (node.lastStart() + 1 == next.firstStart() && node.lastEnd() + 1 == next.firstEnd()) {
                    if (candidate != null) {
                        // found a second candidate
                        return null;
                    }
                    candidate = next;
                }
            }
        }
        return candidate;
    }

    public T getUniquePredecessor(T node) {
        T candidate = null;
        for (long kmer : KmerEncodingHelper.prevStates(k, node.firstKmer())) {
            Object lookup = endKmerLookup.get(kmer);
            if (lookup == null) {
                // nothing to do
            } else if (lookup instanceof Int2ObjectSortedMap) {
                Int2ObjectSortedMap<T> candidates = ((Int2ObjectSortedMap<T>)lookup);
                T prev = candidates.get(node.lastEnd() - 1);
                if (prev != null && prev.lastStart() + 1 == node.firstStart()) {
                    if (candidate != null) {
                        // found a second candidate
                        return null;
                    }
                    candidate = prev;
                }
            } else {
                T prev = (T)lookup;
                if (prev.lastStart() + 1 == node.firstStart() && prev.lastEnd() + 1 == node.firstEnd()) {
                    if (candidate != null) {
                        // found a second candidate
                        return null;
                    }
                    candidate = prev;
                }
            }
        }
        return candidate;
    }
    public void remove(T kn) {
        removeFromStartKmerLookup(kn);
        removeFromEndKmerLookup(kn);
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
            startKmerLookup.remove(kmer);
        }
    }
    public void removeFromEndKmerLookup(T kn) {
        long kmer = kn.lastKmer();
        Object x = endKmerLookup.get(kmer);
        if (x instanceof Int2ObjectSortedMap) {
            Int2ObjectSortedMap<T> positionLookup = (Int2ObjectSortedMap<T>) x;
            T found = positionLookup.remove(kn.lastEnd());
            assert (found != null);
            if (positionLookup.isEmpty()) {
                endKmerLookup.remove(kmer);
            }
        } else {
            endKmerLookup.remove(kmer);
        }
    }
    public void add(T kn) {
        addToStartKmerLookup(kn);
        addOrReplaceEndKmerLookup(kn, kn);
        size++;
    }

    private void addToStartKmerLookup(T kn) {
        long kmer = kn.firstKmer();
        Object x = startKmerLookup.get(kmer);
        Int2ObjectSortedMap<T> positionLookup;
        if (x == null) {
            startKmerLookup.put(kmer, kn);
            return;
        } else if (x instanceof Int2ObjectSortedMap) {
            positionLookup = (Int2ObjectSortedMap<T>) x;
        } else {
            positionLookup = new Int2ObjectRBTreeMap<>();
            startKmerLookup.put(kmer, positionLookup);
            T existing = (T)x;
            positionLookup.put(existing.firstEnd(), existing);
        }
        positionLookup.put(kn.firstEnd(), kn);
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
            if (existing.lastEnd() == key.lastEnd()) {
                // We are doing a replace
                // since we only have 1 record, we can just overwrite it
                endKmerLookup.put(kmer, value);
                return;
            }
            positionLookup = new Int2ObjectRBTreeMap<>();
            endKmerLookup.put(kmer, positionLookup);
            positionLookup.put(existing.lastEnd(), existing);
        }
        positionLookup.put(key.lastEnd(), value);
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

    public void sanityCheck() {
        assert(startKmerLookup.size() == endKmerLookup.size());

        List<KmerNode> start = startKmerLookup.values()
                .stream()
                .flatMap(x -> (x instanceof Int2ObjectSortedMap) ? ((Int2ObjectSortedMap<T>)x).values().stream() : Stream.of((T)x))
                .collect(Collectors.toList());
        List<KmerNode> end = endKmerLookup.values()
                .stream()
                .flatMap(x -> (x instanceof Int2ObjectSortedMap) ? ((Int2ObjectSortedMap<T>)x).values().stream() : Stream.of((T)x))
                .collect(Collectors.toList());
        assert(start.size() == end.size());
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
                assert(((Int2ObjectSortedMap<T>)x).get(kn.lastEnd()) == kn);
            } else {
                assert(x == kn);
            }
        }
    }
}
