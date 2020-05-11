package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.IntegerIntervalSet;
import au.edu.wehi.idsv.util.IntervalUtil;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;


public class KmerPathNodePath extends KmerPathNodeBasePath {
	private ArrayDeque<TraversalNode> nodepath = new ArrayDeque<TraversalNode>();
	private ArrayDeque<KmerPathNode> path = new ArrayDeque<KmerPathNode>();
	private ArrayDeque<Iterator<TraversalNode>> nextPath = new ArrayDeque<Iterator<TraversalNode>>();
	public KmerPathNodePath(KmerPathSubnode node, boolean traverseForward, int maxPathLength) {
		super(node, maxPathLength, traverseForward);
		push(rootNode());
	}
	public int pathLength() {
		return (traversingForward() ? nodepath.getLast() : nodepath.getFirst()).pathLength();
	}
	public int pathWeight() {
		return headNode().pathWeight();
	}
	protected KmerPathNode headPath() {
		return traversingForward() ? path.getLast() : path.getFirst();
	}
	protected TraversalNode headNode() {
		return traversingForward() ? nodepath.getLast() : nodepath.getFirst();
	}
	protected Iterator<TraversalNode> headNext() {
		return traversingForward() ? nextPath.getLast() : nextPath.getFirst();
	}
	/**
	 * Traverse to the next child of this node
	 * @return
	 */
	public boolean dfsNextChild() {
		if (headNext().hasNext()) {
			push(headNext().next());
			return true;
		}
		return false;
	}
	public void push(TraversalNode next) {
		if (traversingForward()) {
			nodepath.addLast(next);
			path.addLast(next.node().node());
			nextPath.addLast(next.iterator());
		} else {
			nodepath.addFirst(next);
			path.addFirst(next.node().node());
			nextPath.addFirst(next.iterator());
		}
	}
	/**
	 * Resets traversal of all child nodes to an untraversed state
	 */
	public void dfsResetChildTraversal() {
		TraversalNode node = headNode();
		popUnchecked();
		push(node);
	}
	/**
	 * Stop any further child traversal of this node and remove it from the path
	 */
	public void pop() {
		if (nodepath.size() == 1) throw new IllegalStateException("Cannot remove root node from traversal path");
		popUnchecked();
	}
	private void popUnchecked() {
		if (traversingForward()) {
			nodepath.removeLast();
			path.removeLast();
			nextPath.removeLast();
		} else {
			nodepath.removeFirst();
			path.removeFirst();
			nextPath.removeFirst();
		}
	}
	public Deque<KmerPathNode> currentPath() {
		return path;
	}
	public IntegerIntervalSet terminalRanges() {
		return headNode().terminalRanges();
	}
	public IntegerIntervalSet terminalLeafRanges() {
		return headNode().terminalLeafAnchorRanges();
	}
	public void greedyTraverse(boolean allowReference, boolean allowNonReference) {
		TraversalNode best;
		do {
			best = null;
			Iterator<TraversalNode> it = headNext();
			while (it.hasNext()) {
				TraversalNode n = it.next();
				boolean isRef = n.node().node().isReference();
				if ((isRef && allowReference) || (!isRef && allowNonReference)) {
					if (best == null || n.node().weight() > best.node().weight()) {
						best = n;
					}
				}
			}
			if (best != null) {
				push(best);
			}
		} while (best != null);
	}
	public void push(Iterator<KmerPathSubnode> pathIt) {
		while (pathIt.hasNext()) {
			KmerPathSubnode sn = pathIt.next();
			push(sn);
		}
	}
	public void push(KmerPathSubnode toNode) {
		Iterator<TraversalNode> it = headNext();
		while (it.hasNext()) {
			TraversalNode tn = it.next();
			if (tn.node().node() == toNode.node()) {
				if (IntervalUtil.overlapsClosed(tn.node().firstStart(), tn.node().firstEnd(), toNode.firstStart(), toNode.firstEnd())) {
					push(tn);
					return;
				}
			}
		}
		throw new IllegalStateException("Illegal traversal");
	}
	public String toString() {
		return headNode().asSubnodes().toString().replace(",", "\n");
	}
}
