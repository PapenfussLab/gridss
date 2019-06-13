package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNodeUtil;

import java.util.*;

/**
 * Backed by ArrayList under the assumption that the number of entries
 * at each position will be small
 */
public class KmerNodeByFirstStartEndKmerReferenceNavigablePartiallyOrderedSet<T extends KmerNode> extends SortedByPositionArrayBackedNavigablePartiallyOrderedSet<T> {
    public KmerNodeByFirstStartEndKmerReferenceNavigablePartiallyOrderedSet(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(T obj) {
        return obj.firstStart();
    }

    @Override
    protected boolean areConsideredEqual(T o1, T o2) {
        return o1.firstEnd() ==  o2.firstEnd()
                && o1.firstKmer() == o2.firstKmer()
                && o1.isReference() == o2.isReference();
    }

    @Override
    public Comparator<? super T> comparator() {
        return KmerNodeUtil.ByFirstStartEndKmerReference;
    }
}
