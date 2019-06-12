package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Queue;
import java.util.stream.Stream;

/**
 * Implements a priority queue keyed only by position.
 *
 * Use this optimised data structure if the order in which records with the same position are returned does not matter
 * and records are only ever removed from the head of the queue.
 */
public abstract class SortedByPositionUnorderedWithinPosition<T> extends SortedByPositionCollection<T, ArrayDeque<T>> implements Queue<T> {
    public SortedByPositionUnorderedWithinPosition(int blockBits) {
        super(blockBits);
    }

    @Override
    protected T peekAtPosition(ArrayDeque<T> ts) {
        return ts.peekFirst();
    }

    @Override
    protected T popAtPosition(ArrayDeque<T> coll) {
        return coll.pollFirst();
    }

    @Override
    protected ArrayDeque<T> createAtPosition() {
        return new ArrayDeque<>();
    }

    @Override
    public boolean offer(T t) {
        return add(t);
    }

    @Override
    public boolean remove(Object x) {
        throw new UnsupportedOperationException("Different different data structure required for efficient remove()");
    }
}
