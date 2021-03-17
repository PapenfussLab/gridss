package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerNode;

public class KmerNodeByLastKmerIntervalLookup<T extends KmerNode> extends KmerIntervalLookup<T> {
    @Override
    protected int getStart(T node) { return node.lastStart(); }

    @Override
    protected int getEnd(T node) {
        return node.lastEnd();
    }

    @Override
    protected long getKmer(T node) {
        return node.lastKmer();
    }
}
