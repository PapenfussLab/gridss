package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerNode;

public class KmerNodeByFirstStartPriorityQueue<T extends KmerNode> extends SortedByPositionUnorderedWithinPosition<T> {
    public KmerNodeByFirstStartPriorityQueue(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(T obj) {
        return obj.firstStart();
    }
}
