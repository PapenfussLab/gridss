package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.TraversalNode;
import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectSortedMap;

import java.util.*;
import java.util.stream.Stream;

public class TraversalNodeByScoreDescPathFirstIdentity implements SortedSet<TraversalNode> {
    private final Int2ObjectSortedMap<Int2ObjectSortedMap<TreeSet<TraversalNode>>> scorePathFirstLookup = new Int2ObjectRBTreeMap<>();
    private int size = 0;
    @Override
    public Comparator<? super TraversalNode> comparator() {
        return TraversalNode.ByScoreDescPathFirstEndSubnode;
    }

    @Override
    public SortedSet<TraversalNode> subSet(TraversalNode fromElement, TraversalNode toElement) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> headSet(TraversalNode toElement) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> tailSet(TraversalNode fromElement) {
        throw new UnsupportedOperationException();
    }

    @Override
    public TraversalNode first() {
        if (size == 0) throw new NoSuchElementException();
        // TODO: compare performance: iterator creation vs multiple tree traversals
        //TreeSet<TraversalNode> lookup = scorePathFirstLookup.values().iterator().next().values().iterator().next();
        Int2ObjectSortedMap<TreeSet<TraversalNode>> pathLookup = scorePathFirstLookup.get(scorePathFirstLookup.firstIntKey());
        TreeSet<TraversalNode> lookup = pathLookup.get(pathLookup.firstIntKey());
        return lookup.first();
    }

    @Override
    public TraversalNode last() {
        throw new UnsupportedOperationException();
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size == 0;
    }

    @Override
    public boolean contains(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<TraversalNode> iterator() {
        return stream().iterator();
    }

    public Stream<TraversalNode> stream() {
        return scorePathFirstLookup.values().stream()
                .flatMap(x -> x.values().stream()
                        .flatMap(y -> y.stream()));
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException();
    }
    private static volatile int maxSize = 0;
    @Override
    public boolean add(TraversalNode tn) {
        Int2ObjectSortedMap<TreeSet<TraversalNode>> pathLookup = scorePathFirstLookup.get(-tn.score);
        if (pathLookup == null) {
            pathLookup = new Int2ObjectRBTreeMap<>();
            scorePathFirstLookup.put(-tn.score, pathLookup);
        }
        TreeSet<TraversalNode> lookup = pathLookup.get(tn.pathFirstStart());
        if (lookup == null) {
            lookup = new TreeSet<>(TraversalNode.BySubnode);
            pathLookup.put(tn.pathFirstStart(), lookup);
        }
        if (lookup.add(tn)) {
            size++;
            return true;
        }
        return false;
    }

    @Override
    public boolean remove(Object o) {
        if (!(o instanceof TraversalNode)) return false;
        return remove((TraversalNode)o);
    }

    public boolean remove(TraversalNode tn) {
        Int2ObjectSortedMap<TreeSet<TraversalNode>> pathLookup = scorePathFirstLookup.get(-tn.score);
        if (pathLookup == null) return false;
        TreeSet<TraversalNode> lookup = pathLookup.get(tn.pathFirstStart());
        if (lookup == null) return false;
        if (lookup.remove(tn)) {
            size--;
            if (lookup.size() == 0) {
                pathLookup.remove(tn.pathFirstStart());
                if (pathLookup.size() == 0) {
                    scorePathFirstLookup.remove(-tn.score);
                }
            }
            return true;
        }
        return false;
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(Collection<? extends TraversalNode> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        boolean changed = false;
        for (Object o : c) {
            changed |= remove(o);
        }
        return changed;
    }

    @Override
    public void clear() {
        size = 0;
        scorePathFirstLookup.clear();
    }
}
