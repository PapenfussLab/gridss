package au.edu.wehi.idsv.debruijn.positional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import au.edu.wehi.idsv.util.IntervalUtil;


/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * Paths consist of a sequence of non-reference nodes with
 * optimal starting and ending reference node anchors.
 * 
 * At each node, the best path for each position interval
 * from each parent node is memoized.
 * 
 * Key to the memoization, is the property that for each
 * unmemoized path, such a path must have started at a
 * starting node with no valid predecessor. To fully
 * memoize the graph a start-coordinate frontier of
 * partial paths can be constructed. Initially populated
 * with all starting nodes, The next nodes from the 
 * frontier head is are traversed and memoized, with
 * each child path added back into the frontier if it
 * is the optimal path to that node.
 * 
 * Since reference nodes can only form the first or last
 * node on a path, they are treated as a special case.
 * - Reference nodes always start a new path across their
 *  entire length.
 * - Path consisting only of reference nodes are not
 *  included in the scoring
 * - Reference nodes always terminate any path leading
 *  to them. 
 * 
 * @author Daniel Cameron
 *
 */
public class MemoizedContigCaller extends ContigCaller {
	/**
	 * Path scores in order of descending score
	 */
	private final SortedSet<TraversalNode> contigByScore = new TreeSet<>(TraversalNode.ByScoreDescPathFirstStartArbitrary);
	private final SortedSet<TraversalNode> frontierByPathStart = new TreeSet<>(TraversalNode.ByPathFirstStartArbitrary);
	private final MemoizedContigTraverse frontier = new MemoizedContigTraverse();
	private int maxVisitedEndPosition = Integer.MIN_VALUE;
	private class MemoizedContigTraverse extends MemoizedTraverse {
		@Override
		protected void onMemoizeAdd(TraversalNode tn) {
			assert(tn.score != ANCHORED_SCORE * 2); // shouldn't ever have reference-reference paths
			if (tn.score == ANCHORED_SCORE) {
				// don't track reference-only paths
				assert(tn.node.isReference());
				assert(tn.parent == null);
				return;
			}
			contigByScore.add(tn);
		}
		@Override
		protected void onMemoizeRemove(TraversalNode tn) {
			contigByScore.remove(tn);
		}
		@Override
		protected void onFrontierAdd(TraversalNode tn) {
			frontierByPathStart.add(tn);
		}
		@Override
		protected void onFrontierRemove(TraversalNode tn) {
			frontierByPathStart.remove(tn);
		}
	}
	public MemoizedContigCaller(
			Iterator<KmerPathNode> it,
			int maxEvidenceWidth) {
		super(it, maxEvidenceWidth);
	}
	/**
	 * Adds a new node to the graph.
	 * 
	 * Parent nodes predecessors need to be added to frontier
	 * so paths to the added node can be calculated
	 * 
	 * @param node node to removed
	 */
	@Override
	public void add(KmerPathNode node) {
		assert(frontier.memoized(node).isEmpty());
		TraversalNode tn = new TraversalNode(new KmerPathSubnode(node), node.isReference() ? ANCHORED_SCORE : 0);
		frontier.memoize(tn);
		if (node.firstStart() > maxVisitedEndPosition + 1) {
			// don't need to flag parents for revisitation
			// since we have never visited any of them
		} else {
			// flag predecessors for recalculation so this node will be visited
			for (KmerPathNode prev : node.prev()) {
				for (TraversalNode prevtn : frontier.memoized(prev)) {
					if (IntervalUtil.overlapsClosed(prevtn.node.lastStart() + 1, prevtn.node.lastEnd() + 1, node.firstStart(), node.firstEnd())) {
						frontier.addFrontier(prevtn);
					}
				}
			}
		}
	}
	/**
	 * Advance the frontier as far as possible
	 */
	private void advanceFrontier() {
		// Can only advance frontier if all possible successors are guaranteed to be loaded
		while (!frontier.isEmptyFrontier() && frontier.peekFrontier().node.lastEnd() < nextPosition() - 1) {
			visit(frontier.pollFrontier());
		}
	}
	/**
	 * Visits the given frontier node
	 * @param node node to visit
	 */
	private void visit(TraversalNode node) {
		assert(node.node.lastEnd() + 1 < nextPosition()); // successors must be fully defined
		if (node.node.isReference()) {
			// Reset reference traversal to a starting path
			node = new TraversalNode(new KmerPathSubnode(node.node.node()), ANCHORED_SCORE);
		}
		for (KmerPathSubnode sn : node.node.next()) {
			if (sn.node().isReference() && node.node.isReference()) {
				// drop reference - reference transitions as we're only
				// traversing non-reference nodes
			} else {
				TraversalNode tn = new TraversalNode(node, sn, sn.isReference() ? ANCHORED_SCORE : 0);
				frontier.memoize(tn);
			}
		}
		maxVisitedEndPosition = Math.max(maxVisitedEndPosition, node.node.lastEnd());
		// We don't need to explicitly track terminal paths.
		// Since node weights are strictly positive, any successor nodes
		// will result in a path with higher score than the interval
		// that terminates here thus a TraversalNode with any successors
		// will never be the highest scoring path
	}
	/**
	 * Removes a node from the graph.
	 * 
	 * @param node node to removed
	 */
	@Override
	public void remove(KmerPathNode node) {
		frontier.remove(node);
		// flag the successors as potential starting nodes
		for (KmerPathNode child : node.next()) {
			TraversalNode tn = new TraversalNode(new KmerPathSubnode(child), child.isReference() ? ANCHORED_SCORE : 0);
			frontier.memoize(tn);
		}
	}
	/**
	 * Determines whether the current best contig is
	 * the globally best contig containing the evidence.
	 * 
	 *  ---- best contig
	 *      ****** overlapping read (maxEvidenceWidth)
	 *            -------------------- non-terminal path
	 *                                | nextPosition()
	 *            ^- earliest non-terminal path
	 *
	 * Can only call best contig when at least maxEvidenceWidth
	 * bases exist between the end of the contig and the start of 
	 * the closest incomplete contig.
	 */
	private boolean canCallBestContig() {
		if (contigByScore.isEmpty()) return false;
		int unprocessedPosition = nextPosition();
		if (!frontierByPathStart.isEmpty()) {
			int frontierPathFirstStart = frontierByPathStart.first().pathFirstStart();
			unprocessedPosition = Math.min(unprocessedPosition, frontierPathFirstStart);
		}
		int bestContigLastEnd = contigByScore.first().node.lastEnd(); 
		return bestContigLastEnd < unprocessedPosition - maxEvidenceWidth - 1;
	}
	@Override
	public ArrayDeque<KmerPathSubnode> bestContig() {
		advanceFrontier();
		while (underlying.hasNext() && !canCallBestContig()) {
			// by loading all nodes within maxEvidenceDistance, we guarantee
			// that all final nodes of incomplete memoized paths
			// can be fully memoized.
			int loadBefore = nextPosition() + 2 * maxEvidenceWidth;
			while (underlying.hasNext() && nextPosition() <= loadBefore) {
				KmerPathNode n = underlying.next();
				add(n);
			}
			advanceFrontier();
		}
		if (contigByScore.isEmpty()) return null;
		TraversalNode tn = contigByScore.first();
		ArrayDeque<KmerPathSubnode> contig = tn.toSubnodeNextPath();
		if (contig.peekFirst().isReference()) {
			contig.pollFirst();
		}
		if (contig.peekLast().isReference()) {
			contig.pollLast();
		}
		assert(!contig.isEmpty());
		assert(contig.stream().allMatch(sn -> !sn.isReference()));
		return contig;
	}
	public int tracking_memoizedNodeCount() {
		return frontier.tracking_memoizedNodeCount();
	}
	public int tracking_frontierSize() {
		return frontier.tracking_frontierSize();
	}
	@Override
	public void exportState(File file) throws IOException {
		frontier.export(file);
	}
}
