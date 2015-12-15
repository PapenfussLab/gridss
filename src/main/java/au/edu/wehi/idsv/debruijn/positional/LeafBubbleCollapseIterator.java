package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

import com.google.common.collect.Range;

/**
 * Collapses leaves and bubbles
 * 
 * 
 * @author Daniel Cameron
 *
 */
public class LeafBubbleCollapseIterator extends CollapseIterator {
	/**
	 * Additional buffer size allowance for chaining of path collapse.
	 * This allows for collapsing of underlying bubbles that have
	 * bubbles branching off them (the underlying bubble is only
	 * considered a bubble after the first collapse)
	 *         B            
	 *        / \
	 *       A - A      Example: can't collapse A until B is collapsed
	 *      /     \
	 *  * - * - * - * - *
	 */
	private static final int RECOLLAPSE_MARGIN = 1;
	public LeafBubbleCollapseIterator(
			Iterator<KmerPathNode> it,
			int k,
			int maxPathCollapseLength,
			int maxBasesMismatch) {
		super(it, k, maxPathCollapseLength, maxBasesMismatch, 0, RECOLLAPSE_MARGIN * maxPathCollapseLength);
	}
	@Override
	protected boolean collapse(KmerPathNode node, int maxCollapseLength) {
		if (processForward(node, maxCollapseLength)) return true;
		if (processBackward(node, maxCollapseLength)) return true;
		return false;
	}
	private boolean processBackward(KmerPathNode node, int maxCollapseLength) {
		for (KmerPathSubnode startCandidate : new KmerPathSubnode(node).subnodesOfDegree(KmerPathSubnode.NOT_MULTIPLE_EDGES, KmerPathSubnode.SINGLE_EDGE)) {
			KmerPathSubnode rootCandidate = startCandidate.next().get(0);
			if (node != rootCandidate.node()) { // don't collapse self loops
				for (Range<Integer> r : rootCandidate.prevPathRangesOfDegree(KmerPathSubnode.MULTIPLE_EDGES).asRanges()) {
					KmerPathSubnode leafStart = new KmerPathSubnode(startCandidate.node(), r.lowerEndpoint() - startCandidate.length(), r.upperEndpoint() - startCandidate.length());
					Set<KmerPathNode> visited = Collections.newSetFromMap(new IdentityHashMap<KmerPathNode, Boolean>());
					visited.add(leafStart.node());
					visited.add(rootCandidate.node());
					if (backwardLeafTraverse(visited, new TraversalNode(new TraversalNode(rootCandidate, 0), leafStart), maxCollapseLength)) return true;
				}
			}
		}
		return false;
	}
	private boolean processForward(KmerPathNode node, int maxCollapseLength) {
		for (KmerPathSubnode startCandidate : new KmerPathSubnode(node).subnodesOfDegree(KmerPathSubnode.SINGLE_EDGE, KmerPathSubnode.NOT_MULTIPLE_EDGES)) {
			KmerPathSubnode rootCandidate = startCandidate.prev().get(0);
			if (node != rootCandidate.node()) { // don't collapse self loops
				for (Range<Integer> r : rootCandidate.nextPathRangesOfDegree(KmerPathSubnode.MULTIPLE_EDGES).asRanges()) {
					KmerPathSubnode leafStart = new KmerPathSubnode(startCandidate.node(), r.lowerEndpoint() + rootCandidate.length(), r.upperEndpoint() + rootCandidate.length());
					Set<KmerPathNode> visited = Collections.newSetFromMap(new IdentityHashMap<KmerPathNode, Boolean>());
					visited.add(leafStart.node());
					visited.add(rootCandidate.node());
					if (forwardLeafTraverse(visited, new TraversalNode(new TraversalNode(rootCandidate, 0), leafStart), maxCollapseLength)) return true;
				}
			}
		}
		return false;
	}
	/**
	 * Finds all backward leaf paths and attempts to merge
	 * 
	 * @param tn leaf path
	 * @return true if a path could be merged, false otherwise
	 */
	private boolean backwardLeafTraverse(Set<KmerPathNode> visited, TraversalNode tn, int maxCollapseLength) {
		KmerPathSubnode node = tn.node;
		for (Range<Integer> range : node.prevPathRangesOfDegree(KmerPathSubnode.NO_EDGES).asRanges()) {
			// Terminal leaf
			TraversalNode terminalNode = new TraversalNode(tn, range.lowerEndpoint(), range.upperEndpoint());
			if (memoizedCollapse(visited, terminalNode, false, null)) return true;
		}
		for (Range<Integer> range : node.prevPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE).asRanges()) {
			KmerPathSubnode sn = new KmerPathSubnode(node.node(), range.lowerEndpoint(), range.upperEndpoint());
			KmerPathSubnode adjNode = sn.prev().get(0);
			if (tn.pathLength + adjNode.length() <= maxCollapseLength && !visited.contains(adjNode.node())) {
				for (Range<Integer> adjRange : adjNode.nextPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE).asRanges()) {
					KmerPathSubnode adjsn = new KmerPathSubnode(adjNode.node(), adjRange.lowerEndpoint(), adjRange.upperEndpoint());
					TraversalNode adjTraveral = new TraversalNode(tn, adjsn);
					visited.add(adjNode.node());
					if (backwardLeafTraverse(visited, adjTraveral, maxCollapseLength)) return true;
					visited.remove(adjNode.node());
				}
				for (Range<Integer> adjRange : adjNode.nextPathRangesOfDegree(KmerPathSubnode.MULTIPLE_EDGES).asRanges()) {
					// End of bubble
					KmerPathSubnode adjsn = new KmerPathSubnode(adjNode.node(), adjRange.lowerEndpoint(), adjRange.upperEndpoint());
					TraversalNode adjTraveral = new TraversalNode(tn, adjsn);
					visited.add(adjNode.node());
					if (memoizedCollapse(visited, adjTraveral, false, adjTraveral.node.node())) return true;
					visited.remove(adjNode.node());
				}
			}
		}
		return false;
	}
	private boolean forwardLeafTraverse(Set<KmerPathNode> visited, TraversalNode tn, int maxCollapseLength) {
		KmerPathSubnode node = tn.node;
		for (Range<Integer> range : node.nextPathRangesOfDegree(KmerPathSubnode.NO_EDGES).asRanges()) {
			// Terminal leaf
			TraversalNode terminalNode = new TraversalNode(tn, range.lowerEndpoint(), range.upperEndpoint());
			if (memoizedCollapse(visited, terminalNode, true, null)) return true;
		}
		for (Range<Integer> range : node.nextPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE).asRanges()) {
			KmerPathSubnode sn = new KmerPathSubnode(node.node(), range.lowerEndpoint(), range.upperEndpoint());
			KmerPathSubnode adjNode = sn.next().get(0);
			if (tn.pathLength + adjNode.length() <= maxCollapseLength && !visited.contains(adjNode.node())) {
				for (Range<Integer> adjRange : adjNode.prevPathRangesOfDegree(KmerPathSubnode.SINGLE_EDGE).asRanges()) {
					KmerPathSubnode adjsn = new KmerPathSubnode(adjNode.node(), adjRange.lowerEndpoint(), adjRange.upperEndpoint());
					TraversalNode adjTraveral = new TraversalNode(tn, adjsn);
					visited.add(adjNode.node());
					if (forwardLeafTraverse(visited, adjTraveral, maxCollapseLength)) return true;
					visited.remove(adjNode.node());
				}
			}
		}
		return false;
	}
	private static class MemoizedPath {
		public final int basesDifferent;
		public final TraversalNode path;
		public MemoizedPath(TraversalNode path, int basesDifferent) {
			this.path = path;
			this.basesDifferent = basesDifferent;
		}
	}
	private void merge(TraversalNode toCollapse, TraversalNode into, boolean traversalForward) {
		if (traversalForward) {
			int commonStart = Math.max(toCollapse.node.lastStart() - toCollapse.pathLength, into.node.lastStart() - into.pathLength);
			int commonEnd = Math.min(toCollapse.node.lastEnd() - toCollapse.pathLength, into.node.lastEnd() - into.pathLength);
			leavesCollapsed++;
			mergeForward(new TraversalNode(toCollapse, commonStart + toCollapse.pathLength - toCollapse.node.length() + 1, commonEnd + toCollapse.pathLength - toCollapse.node.length() + 1),
					new TraversalNode(into, commonStart + into.pathLength - into.node.length() + 1, commonEnd + into.pathLength - into.node.length() + 1));
		} else {
			// take the minimal overlap
			int anchorStart = Math.max(toCollapse.node.firstStart() + toCollapse.pathLength, into.node.firstStart() + into.pathLength);
			int anchorEnd = Math.min(toCollapse.node.firstEnd() + toCollapse.pathLength, into.node.firstEnd() + into.pathLength);
			mergeBackward(new TraversalNode(toCollapse, anchorStart - toCollapse.pathLength, anchorEnd - toCollapse.pathLength),
					new TraversalNode(into, anchorStart - into.pathLength, anchorEnd - into.pathLength));
		}
	}
	private List<TraversalNode> successors(Set<KmerPathNode> refnodes, TraversalNode node, boolean traversalForward) {
		List<TraversalNode> succ = new ArrayList<TraversalNode>(4);
		for (KmerPathSubnode sn : traversalForward ? node.node.next() : node.node.prev()) {
			if (!intersects(refnodes, node, sn)) {
				TraversalNode succtn = new TraversalNode(node, sn);
				succ.add(succtn);
			}
		}
		return succ;
	}
	private static boolean intersects(Set<KmerPathNode> refnodes, TraversalNode node, KmerPathSubnode sn) {
		if (refnodes.contains(sn.node())) return true;
		return node.traversingWouldCauseSelfIntersection(sn.node());
	}
	private int partialSequenceBasesDifferent(LongArrayList toCollapsePathKmers, TraversalNode tn, boolean traversalForward) {
		LongArrayList nodeKmers = tn.node.node().pathKmers();
		int basesDifference;
		if (traversalForward) {
			basesDifference = KmerEncodingHelper.partialSequenceBasesDifferent(k, toCollapsePathKmers, nodeKmers, tn.pathLength - tn.node.length(), true);
		} else {
			basesDifference = KmerEncodingHelper.partialSequenceBasesDifferent(k, toCollapsePathKmers, nodeKmers, toCollapsePathKmers.size() - tn.pathLength, false);
		}
		return basesDifference;
	}
	private static List<MemoizedPath> frontierPop(SortedMap<KmerPathSubnode, List<MemoizedPath>> frontier) {
		KmerPathSubnode key = frontier.firstKey();
		List<MemoizedPath> values = frontier.get(key);
		frontier.remove(key);
		return values;
	}
	private boolean frontierProcess(SortedMap<KmerPathSubnode, List<MemoizedPath>> frontier, MemoizedPath mp, TraversalNode toCollapse, boolean traversalForward, KmerPathNode terminalNode) {
		if (mp.basesDifferent > maxBasesMismatch) return false;
		if (mp.path.pathLength >= toCollapse.pathLength) {
			// terminal node
			if (mp.path.score < toCollapse.score) return false;
			// leaf = collapse
			// branch = must end at same node
			if (terminalNode != null && (mp.path.node.node() != terminalNode || mp.path.pathLength != toCollapse.pathLength)) {
				// bubble must end at the correct node with matching path length
				return false;
			}
			if (terminalNode == null) {
				leavesCollapsed++;
			} else {
				branchesCollapsed++;
			}
			merge(toCollapse, mp.path, traversalForward);
			return true;
		}
		frontierAdd(frontier, mp);
		return false;
	}
	/**
	 * Adds the given node to the frontier
	 * @param frontier frontier
	 * @param mp node to add
	 */
	private static void frontierAdd(SortedMap<KmerPathSubnode, List<MemoizedPath>> frontier, MemoizedPath mp) {
		List<MemoizedPath> mpl = frontier.get(mp.path.node);
		if (mpl == null) {
			mpl = new ArrayList<MemoizedPath>(3);
			frontier.put(mp.path.node, mpl);
		}
		for (int i = 0; i < mpl.size(); i++) {
			MemoizedPath n = mpl.get(i);
			if (n.basesDifferent >= mp.basesDifferent && n.path.score < mp.path.score) {
				// existing path is more different, with less weight = ours is better
				mpl.remove(i);
				i--;
			} else if (n.basesDifferent <= mp.basesDifferent && n.path.score >= mp.path.score) {
				// existing path is less different, with better score = ours is worse
				return;
			}
		}
		mpl.add(mp);
	}
	/**
	 * Use a memoized breadth first search to find a similar path
	 * @param collapseNodes lookup of collapse path
	 * @param toCollapse path to attempt to collapse
	 * @param traversalForward traversal direction
	 * @param terminalNode node our path must finish on, null if collapsing leaf
	 * @return true if the path was collapsed, false otherwise
	 */
	private boolean memoizedCollapse(Set<KmerPathNode> collapseNodes, TraversalNode toCollapse, boolean traversalForward, KmerPathNode terminalNode) {
		LongArrayList toCollapsePathKmers = new LongArrayList(toCollapse.pathLength);
		for (KmerPathSubnode sn : traversalForward ? toCollapse.toSubnodeNextPath() : toCollapse.toSubnodePrevPath()) {
			toCollapsePathKmers.addAll(sn.node().pathKmers());
		}
		assert(toCollapsePathKmers.size() == toCollapse.pathLength);
		if (terminalNode != null) {
			collapseNodes.remove(terminalNode);
		}
		SortedMap<KmerPathSubnode, List<MemoizedPath>> frontier = new TreeMap<KmerPathSubnode, List<MemoizedPath>>(KmerNodeUtil.ByFirstStartKmer);
		KmerPathSubnode root = traversalForward ? toCollapse.toSubnodeNextPath().getFirst() : toCollapse.toSubnodePrevPath().getLast();
		// set up frontier
		for (TraversalNode tn : successors(collapseNodes, new TraversalNode(root, 0), traversalForward)) {
			MemoizedPath mp = new MemoizedPath(tn, partialSequenceBasesDifferent(toCollapsePathKmers, tn, traversalForward));
			if (frontierProcess(frontier, mp, toCollapse, traversalForward, terminalNode)) return true;
		}
		while (!frontier.isEmpty()) {
			for (MemoizedPath mp : frontierPop(frontier)) {
				for (TraversalNode tn : successors(collapseNodes, mp.path, traversalForward)) {
					int basesDifferent = mp.basesDifferent + partialSequenceBasesDifferent(toCollapsePathKmers, tn, traversalForward);
					MemoizedPath mpnext = new MemoizedPath(tn, basesDifferent);
					if (frontierProcess(frontier, mpnext, toCollapse, traversalForward, terminalNode)) return true;
				}
			}
		}
		if (terminalNode != null) {
			collapseNodes.add(terminalNode);
		}
		return false;
	}
	private void mergeBackward(TraversalNode source, TraversalNode target) {
		int skipSourceCount = Math.max(0, source.pathLength - target.pathLength);
		int skipTargetCount = Math.max(0, target.pathLength - source.pathLength);
		merge(new ArrayList<KmerPathSubnode>(source.toSubnodePrevPath()), new ArrayList<KmerPathSubnode>(target.toSubnodePrevPath()), skipSourceCount, skipTargetCount);
	}
	private void mergeForward(TraversalNode source, TraversalNode target) {
		merge(new ArrayList<KmerPathSubnode>(source.toSubnodeNextPath()), new ArrayList<KmerPathSubnode>(target.toSubnodeNextPath()), 0, 0);
	}
	@Override
	protected boolean reprocessMergedNodes() {
		// still not quite enough to reprocess since leaves branching
		// off the collapsed not are not themselves collapsed since the
		// leaf starting point is considered to be adjacent to the anchor
		// and only the anchor is involved in the merged.
		
		// FIXME: we could get CollapseIterator to add merge neighbours
		// to the unprocessed queue since CollapseIterator guarantees
		// that nodes adjacent to any collapsable path have been loaded
		// into our graph.
		return true;
	}
}