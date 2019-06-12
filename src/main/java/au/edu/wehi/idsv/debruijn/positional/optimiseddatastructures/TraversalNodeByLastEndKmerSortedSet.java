package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.TraversalNode;
import au.edu.wehi.idsv.util.MessageThrottler;
import gridss.ComputeSamTags;
import it.unimi.dsi.fastutil.longs.Long2ObjectRBTreeMap;

import java.util.*;

public class TraversalNodeByLastEndKmerSortedSet extends SortedByPosition<TraversalNode, Long2ObjectRBTreeMap<TraversalNode>> implements SortedSet<TraversalNode> {
    public TraversalNodeByLastEndKmerSortedSet(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(TraversalNode obj) {
        return obj.node.lastEnd();
    }

    @Override
    protected TraversalNode peekAtPosition(Long2ObjectRBTreeMap<TraversalNode> coll) {
        return coll.get(coll.firstLongKey());
    }

    @Override
    protected TraversalNode popAtPosition(Long2ObjectRBTreeMap<TraversalNode> coll) {
        return coll.remove(coll.firstLongKey());
    }

    @Override
    protected Long2ObjectRBTreeMap<TraversalNode> createAtPosition() {
        return new Long2ObjectRBTreeMap<>();
    }

    @Override
    protected boolean addAtPosition(Long2ObjectRBTreeMap<TraversalNode> existing, TraversalNode toAdd) {
        existing.put(toAdd.node.lastKmer(), toAdd);
        return true;
    }

    @Override
    protected boolean removeAtPosition(Long2ObjectRBTreeMap<TraversalNode> coll, TraversalNode obj) {
        return coll.remove(obj.node.lastKmer()) != null;
    }

    @Override
    protected boolean positionIsEmpty(Long2ObjectRBTreeMap<TraversalNode> coll) {
        return coll.isEmpty();
    }

    @Override
    protected boolean containsAtPosition(Long2ObjectRBTreeMap<TraversalNode> coll, TraversalNode obj) {
        return coll.get(obj.node.lastKmer()) != null;
    }

    @Override
    public Comparator<? super TraversalNode> comparator()   {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> subSet(TraversalNode fromElement, TraversalNode toElement)  { throw new UnsupportedOperationException(); }

    @Override
    public SortedSet<TraversalNode> headSet(TraversalNode toElement)  {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> tailSet(TraversalNode fromElement)   {
        throw new UnsupportedOperationException();
    }

    @Override
    public TraversalNode last()  {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<TraversalNode> iterator()  {
        throw new UnsupportedOperationException();
    }

    @Override
    public Spliterator<TraversalNode> spliterator()    {
        throw new UnsupportedOperationException();
    }
}