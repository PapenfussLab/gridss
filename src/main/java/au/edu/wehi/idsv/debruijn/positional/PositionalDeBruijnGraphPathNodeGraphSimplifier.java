package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

public class PositionalDeBruijnGraphPathNodeGraphSimplifier implements Iterator<KmerPathNode> {
	private final PeekingIterator<KmerPathNode> underlying;
	private final int k;
	private final int maxCollapseLength;
	private final int maxBaseMismatchForCollapse;
	private final boolean collapseBubblesOnly;
	private final Queue<KmerPathNode> outputBuffer = new ArrayDeque<KmerPathNode>();
	private final Set<KmerPathNode> active = new HashSet<KmerPathNode>();
	private final PriorityQueue<KmerPathNode> activeBranchStart = new PriorityQueue<KmerPathNode>(1024, KmerPathNode.ByEndPosition);
	private final PriorityQueue<KmerPathNode> activeLeafStartAnchored = new PriorityQueue<KmerPathNode>(1024, KmerPathNode.ByEndPosition);
	private final PriorityQueue<KmerPathNode> activeLeafEndAnchored = new PriorityQueue<KmerPathNode>(1024, KmerPathNode.ByFirstKmerStartPosition);
	private int currentPosition = Integer.MIN_VALUE;
	public PositionalDeBruijnGraphPathNodeGraphSimplifier(Iterator<KmerPathNode> it, int k, int maxCollapseLength, boolean collapseBubblesOnly, int maxBaseMismatchForCollapse) {
		this.underlying = Iterators.peekingIterator(it);
		this.k = k;
		this.maxCollapseLength = maxCollapseLength;
		this.maxBaseMismatchForCollapse = maxBaseMismatchForCollapse;
		this.collapseBubblesOnly = collapseBubblesOnly;
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !outputBuffer.isEmpty();
	}
	@Override
	public KmerPathNode next() {
		ensureBuffer();
		return outputBuffer.poll();
	}
	private void ensureBuffer() {
		if (!outputBuffer.isEmpty()) return;
		// load nodes into graph until we are guaranteed
		// to be able to traverse all branches/leaves
		// rooted at the first node
		if (underlying.hasNext()) {
			currentPosition = underlying.peek().startPosition(0);
			addToGraph();
		} else {
			currentPosition = Integer.MAX_VALUE;
		}
	}
	private void addToGraph() {
		while (underlying.hasNext() && underlying.peek().startPosition(0) <= currentPosition) {
			KmerPathNode pn = underlying.next();
			active.add(pn);
			if (pn.nextList.size() == 0 && pn.prevList.size() == 1) {
				activeLeafStartAnchored.add(pn);
			} else if (pn.prevList.size() == 0 && pn.nextList.size() == 1) {
				activeLeafEndAnchored.add(pn);
			} else if (pn.nextList.size() > 1) {
				activeBranchStart.add(pn);
			}
		}
	}
	private void process() {
		// all paths starting at or before currentPosition are in our graph
		collapseActiveBranches();
		
		//collapseActiveBranchesStartEndingBefore(startPosition - maxCollapseLength - 1);
		//collapseActiveStartAnchoredLeavesEndingBefore(startPosition - maxCollapseLength - 1);
	}
	private void collapseActiveBranches() {
		//while (!activeBranchStart.isEmpty() && activeBranchStart.peek().endPosition() < 
	}
}