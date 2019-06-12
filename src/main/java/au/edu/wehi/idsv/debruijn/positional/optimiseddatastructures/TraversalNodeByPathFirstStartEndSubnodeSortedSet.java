package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.TraversalNode;
import au.edu.wehi.idsv.util.MessageThrottler;
import htsjdk.samtools.util.Log;

import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

public class TraversalNodeByPathFirstStartEndSubnodeSortedSet extends SortedByPosition<TraversalNode, SortedSet<TraversalNode>> implements SortedSet<TraversalNode> {
    private static final Log log = Log.getInstance(TraversalNodeByPathFirstStartEndSubnodeSortedSet.class);


    public TraversalNodeByPathFirstStartEndSubnodeSortedSet(int blockBits) {
        super(blockBits);
    }

    @Override
    protected int getPosition(TraversalNode obj) {
        return obj.pathFirstStart();
    }

    @Override
    protected TraversalNode peekAtPosition(SortedSet<TraversalNode> traversalNodes) {
        return traversalNodes.first();
    }

    @Override
    protected TraversalNode popAtPosition(SortedSet<TraversalNode> traversalNodes) {
        // TODO: check whether this is actually more efficient than first()/remove()
        Iterator<TraversalNode> it = traversalNodes.iterator();
        TraversalNode tn = it.next();
        it.remove();
        return tn;
    }

    @Override
    protected SortedSet<TraversalNode> createAtPosition() {
        return new TreeSet<>(TraversalNode.ByScoreEndSubnode);
    }

    @Override
    protected boolean addAtPosition(SortedSet<TraversalNode> existing, TraversalNode toAdd) {
        return existing.add(toAdd);
    }

    @Override
    protected boolean removeAtPosition(SortedSet<TraversalNode> traversalNodes, TraversalNode obj) {
        return traversalNodes.remove(obj);
    }

    @Override
    protected boolean positionIsEmpty(SortedSet<TraversalNode> traversalNodes) {
        return traversalNodes.isEmpty();
    }

    @Override
    protected boolean containsAtPosition(SortedSet<TraversalNode> traversalNodes, TraversalNode obj) {
        return traversalNodes.contains(obj);
    }

    @Override
    public Comparator<? super TraversalNode> comparator()   {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> subSet(TraversalNode fromElement, TraversalNode toElement)   {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> headSet(TraversalNode toElement)   {
        throw new UnsupportedOperationException();
    }

    @Override
    public SortedSet<TraversalNode> tailSet(TraversalNode fromElement)    {
        throw new UnsupportedOperationException();
    }

    @Override
    public TraversalNode last()   {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<TraversalNode> iterator() {
        if (!MessageThrottler.Current.shouldSupress(log, "performance:TraversalNodeByLastEndKmerSortedSet")) {
            log.warn("TraversalNodeByLastEndKmerSortedSet.iterator() call. This is inefficient and should be no be called in production code.");
        }
        return getCollectionsInOrder().stream().flatMap(ss -> ss.stream()).iterator();
    }
}
