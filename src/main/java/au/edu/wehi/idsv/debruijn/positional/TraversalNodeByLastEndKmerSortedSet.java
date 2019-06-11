package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntRBTreeSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import it.unimi.dsi.fastutil.longs.Long2ObjectSortedMap;

import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedSet;
/*
public class TraversalNodeByLastEndKmerSortedSet implements SortedSet<TraversalNode> {
    private Int2ObjectMap<Long2ObjectSortedMap<TraversalNode>> positions = new Int2ObjectOpenHashMap<>();
    private IntSortedSet occupiedPositions = new IntRBTreeSet();
    private int size = 0;
    @Override
    public Comparator<? super TraversalNode> comparator() {
        return TraversalNode.ByLastEndKmer;
    }

    @Override
    public SortedSet<TraversalNode> subSet(TraversalNode fromElement, TraversalNode toElement) {
        throw new UnsupportedOperationException("NYI");
    }

    @Override
    public SortedSet<TraversalNode> headSet(TraversalNode toElement) {
        throw new UnsupportedOperationException("NYI");
    }

    @Override
    public SortedSet<TraversalNode> tailSet(TraversalNode fromElement) {
        throw new UnsupportedOperationException("NYI");
    }

    @Override
    public TraversalNode first() {
        int firstPosition = occupiedPositions.firstInt();
        Long2ObjectSortedMap<TraversalNode> map = positions.get(firstPosition);
        return map.get(map.firstLongKey());
    }

    @Override
    public TraversalNode last() {
        throw new UnsupportedOperationException("NYI");
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
        throw new UnsupportedOperationException("NYI");
    }

    @Override
    public Iterator<TraversalNode> iterator() {
        throw new UnsupportedOperationException("NYI");
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException("NYI");
    }

    @Override
    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException("NYI");
    }

    @Override
    public boolean add(TraversalNode traversalNode) {
        size++;
    }

    @Override
    public boolean remove(Object o) {
        if (o instanceof TraversalNode) {
            remove((TraversalNode)o);
        } else {
            return false;
        }
    }
    public boolean remove(TraversalNode tn) {
        int position = tn.node.lastEnd();
        Long2ObjectSortedMap<TraversalNode> kmers = positions.get(position);
        long kmer = tn.node.firstKmer();
        if (o instanceof TraversalNode) {
            size--;
        } else {
            return false;
        }
    }

    @Override
    public boolean containsAll(Collection<?> c) {
    }

    @Override
    public boolean addAll(Collection<? extends TraversalNode> c) {
    }

    @Override
    public boolean retainAll(Collection<?> c) {
    }

    @Override
    public boolean removeAll(Collection<?> c) {
    }

    @Override
    public void clear() {
        size = 0;
        positions.clear();
        occupiedPositions.clear();
    }
}
*/