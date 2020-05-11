package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.SortedSet;

/**
 * Backed by ArrayList under the assumption that the number of entries at each position will be small
 */
public abstract class SortedByPositionArrayBackedNavigablePartiallyOrderedSet<T> extends SortedByPosition<T, ArrayList<T>> implements NavigableSet<T> {
    public SortedByPositionArrayBackedNavigablePartiallyOrderedSet(int blockBits) {
        super(blockBits);
    }

    @Override
    protected T peekAtPosition(ArrayList<T> ts) {
        return ts.get(ts.size() - 1);
    }

    @Override
    protected T popAtPosition(ArrayList<T> ts) {
        return ts.remove(ts.size() - 1);
    }

    @Override
    protected ArrayList<T> createAtPosition() {
        return new ArrayList<>(4);
    }

    @Override
    protected boolean addAtPosition(ArrayList<T> existing, T obj) {
        if (containsAtPosition(existing, obj)) return false;
        return existing.add(obj);
    }

    @Override
    protected boolean removeAtPosition(ArrayList<T> existing, T obj) {
        for (int i = 0; i < existing.size(); i++) {
            if (areConsideredEqual(existing.get(i), obj)) {
                existing.remove(i);
                return true;
            }
        }
        return false;
    }

    @Override
    protected boolean positionIsEmpty(ArrayList<T> ts) {
        return ts.isEmpty();
    }

    @Override
    protected int positionSize(ArrayList<T> ts) {
        return ts.size();
    }

    @Override
    protected boolean containsAtPosition(ArrayList<T> existing, T obj) {
        for (int i = 0; i < existing.size(); i++) {
            if (areConsideredEqual(existing.get(i), obj)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Determines if two objects at the given position are to be considered equal
     * @return true if they are the same, false otherwise
     */
    protected abstract boolean areConsideredEqual(T o1, T o2);

    @Override
    public T pollFirst() {
        return poll();
    }

    @Override
    public T lower(T t) {
        throw new UnsupportedOperationException();
    }

    @Override
    public T floor(T t) {
        throw new UnsupportedOperationException();
    }

    @Override
    public T ceiling(T t) {
        throw new UnsupportedOperationException();
    }

    @Override
    public T higher(T t) {
        throw new UnsupportedOperationException();
    }

    @Override
    public T pollLast() {
        throw new UnsupportedOperationException();
    }

    @Override
    public NavigableSet<T> descendingSet() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<T> descendingIterator() {
        throw new UnsupportedOperationException();
    }

    @Override
    public NavigableSet<T> subSet(T t, boolean b, T e1, boolean b1) {
        throw new UnsupportedOperationException();
    }

    @Override
    public NavigableSet<T> headSet(T t, boolean b) {
        throw new UnsupportedOperationException();
    }

    @Override
    public NavigableSet<T> tailSet(T t, boolean b) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<T> subSet(T t, T e1) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<T> headSet(T t) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<T> tailSet(T t) {
        throw new UnsupportedOperationException();
    }

    @Override
    public T last() {
        throw new UnsupportedOperationException();
    }
}
