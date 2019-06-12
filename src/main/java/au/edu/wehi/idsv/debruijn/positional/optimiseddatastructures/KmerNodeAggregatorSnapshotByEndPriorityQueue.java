package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.AggregateNodeIterator;

public class KmerNodeAggregatorSnapshotByEndPriorityQueue extends SortedByPositionUnorderedWithinPosition<AggregateNodeIterator.KmerNodeAggregator.KmerNodeAggregatorSnapshot> {
    public KmerNodeAggregatorSnapshotByEndPriorityQueue(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(AggregateNodeIterator.KmerNodeAggregator.KmerNodeAggregatorSnapshot obj) {
        return obj.snapshotEnd;
    }
}