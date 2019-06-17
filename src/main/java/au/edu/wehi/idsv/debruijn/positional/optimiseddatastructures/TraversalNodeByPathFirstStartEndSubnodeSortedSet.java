package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.TraversalNode;

import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

public class TraversalNodeByPathFirstStartEndSubnodeSortedSet extends SortedByPositionCollection<TraversalNode, SortedSet<TraversalNode>> implements SortedSet<TraversalNode> {
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
}
