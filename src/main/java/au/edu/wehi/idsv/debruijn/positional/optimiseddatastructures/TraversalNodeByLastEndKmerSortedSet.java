package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.TraversalNode;
import it.unimi.dsi.fastutil.longs.Long2ObjectRBTreeMap;

import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedSet;

public class TraversalNodeByLastEndKmerSortedSet extends SortedCollectionByPosition<TraversalNode, Long2ObjectRBTreeMap<TraversalNode>> implements SortedSet<TraversalNode> {
    public TraversalNodeByLastEndKmerSortedSet(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(TraversalNode obj) {
        return obj.node.lastEnd();
    }

    @Override
    protected TraversalNode peekAtPosition(Long2ObjectRBTreeMap<TraversalNode> coll) {
        return coll.get(coll.firstLongKey());
    }

    @Override
    protected TraversalNode popAtPosition(Long2ObjectRBTreeMap<TraversalNode> coll) {
        return coll.remove(coll.firstLongKey());
    }

    @Override
    protected Long2ObjectRBTreeMap<TraversalNode> createAtPosition() {
        return new Long2ObjectRBTreeMap<>();
    }

    @Override
    protected boolean addAtPosition(Long2ObjectRBTreeMap<TraversalNode> existing, TraversalNode toAdd) {
        existing.put(toAdd.node.lastKmer(), toAdd);
        return true;
    }

    @Override
    protected boolean removeAtPosition(Long2ObjectRBTreeMap<TraversalNode> coll, TraversalNode obj) {
        return coll.remove(obj.node.lastKmer()) != null;
    }

    @Override
    protected boolean dataIsEmpty(Long2ObjectRBTreeMap<TraversalNode> coll) {
        return coll.isEmpty();
    }

    @Override
    public Comparator<? super TraversalNode> comparator() {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> subSet(TraversalNode traversalNode, TraversalNode e1) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> headSet(TraversalNode traversalNode) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> tailSet(TraversalNode traversalNode) {
        throw new UnsupportedOperationException();
    }

    @Override
    public TraversalNode first() {
        return peekFirstEntry();
    }

    @Override
    public TraversalNode last() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean contains(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<TraversalNode> iterator() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException();
    }

    @Override
    public <T> T[] toArray(T[] ts) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean containsAll(Collection<?> collection) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(Collection<? extends TraversalNode> collection) {
        boolean changed = false;
        for (TraversalNode n : collection) {
            changed |= add(n);
        }
        return changed;
    }

    @Override
    public boolean retainAll(Collection<?> collection) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean removeAll(Collection<?> collection) {
        boolean changed = false;
        for (Object n : collection) {
            changed |= remove(n);
        }
        return changed;
    }
}