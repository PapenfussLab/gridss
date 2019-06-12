package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerSupportNode;

import java.util.*;

/**
 * Optimised priority queue implementation that efficiently implements only the
 * operations used by SupportNodeIterator
 *
 * Alternate implementation include:
 *  - IntSortedSet lookup backed by Int2ObjectHashMap
 */

public class KmerSupportNodeByFirstStartPriorityQueue extends SortedByPositionCollection<KmerSupportNode, ArrayDeque<KmerSupportNode>> implements Queue<KmerSupportNode> {
    public KmerSupportNodeByFirstStartPriorityQueue(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(KmerSupportNode obj) {
        return obj.firstStart();
    }

    @Override
    protected KmerSupportNode peekAtPosition(ArrayDeque<KmerSupportNode> kmerSupportNodes) {
        return kmerSupportNodes.peekFirst();
    }

    @Override
    protected KmerSupportNode popAtPosition(ArrayDeque<KmerSupportNode> kmerSupportNodes) {
        return kmerSupportNodes.pollFirst();
    }

    @Override
    protected ArrayDeque<KmerSupportNode> createAtPosition() {
        return new ArrayDeque<>();
    }

    @Override
    public boolean offer(KmerSupportNode kmerSupportNode) {
        return add(kmerSupportNode);
    }
}
