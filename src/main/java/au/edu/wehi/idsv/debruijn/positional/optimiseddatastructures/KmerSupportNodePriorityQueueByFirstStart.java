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

public class KmerSupportNodePriorityQueueByFirstStart extends SortedByPosition<KmerSupportNode, ArrayDeque<KmerSupportNode>> implements Queue<KmerSupportNode> {
    public KmerSupportNodePriorityQueueByFirstStart(int blockBits) {
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
    protected boolean addAtPosition(ArrayDeque<KmerSupportNode> existing, KmerSupportNode toAdd) {
        existing.addLast(toAdd);
        return true;
    }

    @Override
    protected boolean removeAtPosition(ArrayDeque<KmerSupportNode> kmerSupportNodes, KmerSupportNode obj) {
        return kmerSupportNodes.remove(obj);
    }

    @Override
    protected boolean positionIsEmpty(ArrayDeque<KmerSupportNode> kmerSupportNodes) {
        return kmerSupportNodes.isEmpty();
    }

    @Override
    protected boolean containsAtPosition(ArrayDeque<KmerSupportNode> kmerSupportNodes, KmerSupportNode obj) {
        return kmerSupportNodes.contains(obj);
    }

    @Override
    public boolean offer(KmerSupportNode kmerSupportNode) {
        return add(kmerSupportNode);
    }

    @Override
    public Iterator<KmerSupportNode> iterator()   {
        throw new UnsupportedOperationException();
    }
}
