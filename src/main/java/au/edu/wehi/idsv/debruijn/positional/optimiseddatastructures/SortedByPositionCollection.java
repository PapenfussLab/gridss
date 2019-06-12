package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import java.util.Collection;
import java.util.stream.Stream;

public abstract class SortedByPositionCollection<T, TColl extends Collection<T>> extends SortedByPosition<T, TColl> implements Collection<T> {
    public SortedByPositionCollection(int blockBits) {
        super(blockBits);
    }

    @Override
    protected T popAtPosition(TColl ts) {
        T t = peekAtPosition(ts);
        remove(t);
        return t;
    }

    @Override
    protected boolean addAtPosition(TColl existing, T toAdd) {
        return existing.add(toAdd);
    }

    @Override
    protected boolean removeAtPosition(TColl ts, T obj) {
        return ts.remove(obj);
    }

    @Override
    protected boolean positionIsEmpty(TColl ts) {
        return ts.isEmpty();
    }

    @Override
    protected int positionSize(TColl ts) {
        return ts.size();
    }

    @Override
    protected boolean containsAtPosition(TColl ts, T obj) {
        return ts.contains(obj);
    }

    @Override
    protected Stream<T> positionStream(TColl coll) {
        return coll.stream();
    }
}
