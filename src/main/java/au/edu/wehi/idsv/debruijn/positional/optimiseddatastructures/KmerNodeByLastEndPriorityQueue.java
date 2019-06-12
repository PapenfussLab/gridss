package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerNode;

public class KmerNodeByLastEndPriorityQueue<T extends KmerNode> extends SortedByPositionUnorderedWithinPosition<T> {
    public KmerNodeByLastEndPriorityQueue(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(T obj) {
        return obj.lastEnd();
    }
}
