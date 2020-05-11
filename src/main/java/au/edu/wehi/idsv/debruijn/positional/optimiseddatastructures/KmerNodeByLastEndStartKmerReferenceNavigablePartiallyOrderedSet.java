package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNodeUtil;

import java.util.Comparator;

/**
 * Backed by ArrayList under the assumption that the number of entries at each position will be small.
 *
 */
public class KmerNodeByLastEndStartKmerReferenceNavigablePartiallyOrderedSet<T extends KmerNode>  extends SortedByPositionArrayBackedNavigablePartiallyOrderedSet<T> {
    public KmerNodeByLastEndStartKmerReferenceNavigablePartiallyOrderedSet(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(T obj) {
        return obj.lastEnd();
    }

    @Override
    protected boolean areConsideredEqual(T o1, T o2) {
        return o1.lastStart() ==  o2.lastStart()
                && o1.lastKmer() == o2.lastKmer()
                && o1.isReference() == o2.isReference();
    }

    @Override
    public Comparator<? super T> comparator() {
        return KmerNodeUtil.ByLastEndStartKmerReference;
    }
}
