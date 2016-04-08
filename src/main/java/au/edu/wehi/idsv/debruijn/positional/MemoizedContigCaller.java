package au.edu.wehi.idsv.debruijn.positional;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.SanityCheckFailureException;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker.MemoizationStats;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.PeekingIterator;


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
	private static final Log log = Log.getInstance(MemoizedContigCaller.class);
	static final boolean ASSERT_ALL_OPERATIONS = false;
	/**
	 * Chunk size (in multiples of maxEvidenceWidth) to load at any one time.
	 * Larger chunks increase memory consumption by reduce the number of
	 * memoization recalculations required.
	 */
	private static final int EVIDENCE_WIDTH_LOADING_MULTIPLE = 2;
	/**
	 * Path scores in order of descending score
	 */
	private final SortedSet<TraversalNode> contigByScore = new TreeSet<>(TraversalNode.ByScoreDescPathFirstEndSubnode);
	private final SortedSet<TraversalNode> frontierByPathStart = new TreeSet<>(TraversalNode.ByPathFirstStartEndSubnode);
	private final MemoizedContigTraverse frontier = new MemoizedContigTraverse();
	private int maxVisitedEndPosition = Integer.MIN_VALUE;
	private class MemoizedContigTraverse extends MemoizedTraverse {
		@Override
		protected void onMemoizeAdd(TraversalNode tn) {
			if (tn.node.isReference() && tn.parent == null) {
				// don't track reference-only paths
				return;
			}
			if (tn.score >= 2 * ANCHORED_SCORE && tn.node.isReference() &&
					tn.parent != null && tn.parent.node.isReference()) {
				// don't track reference-reference paths
				return;
			}
			contigByScore.add(tn);
		}
		@Override
		protected void onMemoizeRemove(TraversalNode tn) {
			contigByScore.remove(tn);
		}
		@Override
		protected void onMemoizeRemove(Collection<TraversalNode> tn) {
			contigByScore.removeAll(tn);
		}
		@Override
		protected void onFrontierAdd(TraversalNode tn) {
			frontierByPathStart.add(tn);
		}
		@Override
		protected void onFrontierRemove(TraversalNode tn) {
			frontierByPathStart.remove(tn);
		}
		@Override
		protected void onFrontierRemove(Collection<TraversalNode> tn) {
			frontierByPathStart.removeAll(tn);
		}
	}
	public MemoizedContigCaller(
			PeekingIterator<KmerPathNode> it,
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
		assert(node.isValid());
		assert(frontier.memoized(node).isEmpty());
		TraversalNode tn = new TraversalNode(new KmerPathSubnode(node), node.isReference() ? ANCHORED_SCORE - node.weight() : 0);
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
		if (Defaults.SANITY_CHECK_MEMOIZATION && ASSERT_ALL_OPERATIONS) {
			sanityCheck();
		}
	}
	/**
	 * Advance the frontier as far as possible. Each frontier advancement must advance
	 * all paths such that no node is in the set of previously visited non-frontier nodes
	 * remaining unvisited after this advancement. If not, the following edge case can
	 * occur:
	 *     
	 *  Frontier: C(ABC), B(DB)
	 *  A - B - C
	 *     /
	 *    D
	 *  
	 *  Remove A:
	 *      B - C
	 *     /
	 *    D
	 *    
	 *   Frontier: C(ABC), B(DB)
	 * Since B does not have any memoization of A, unmemoization stops at B even though
	 * the memoized path to C involves A.
	 * 
	 * Ensuring that the B(DB) frontier is visited will ensure that the orphaned
	 * C(ABC) paths will be removed before invalidation.
	 * 
	 * How can we get into this situation?
	 * For a frontier sorted by lastEnd position, and advancement when the node lastEnd is < unprocessed,
	 * we need node B to have a valid successor C and there to exist a node X such that ... ?
	 * 
	 */
	private void advanceFrontier() {
		// Can only advance frontier if all possible successors are guaranteed to be loaded
		while (!frontier.isEmptyFrontier() && frontier.peekFrontier().node.lastEnd() < nextPosition() - 1) {
			visit(frontier.pollFrontier());
		}
		if (Defaults.SANITY_CHECK_MEMOIZATION && ASSERT_ALL_OPERATIONS) {
			sanityCheck();
			sanityCheckFrontier();
		}
	}
	/**
	 * Visits the given frontier node.
	 * 
	 * @param node node to visit
	 */
	private void visit(TraversalNode node) {
		assert(node.node.lastEnd() + 1 < nextPosition()); // successors must be fully defined
		if (node.node.isReference()) {
			// Reset reference traversal to a starting path
			node = new TraversalNode(new KmerPathSubnode(node.node.node()), ANCHORED_SCORE - node.node.node().weight());
		}
		for (KmerPathSubnode sn : node.node.next()) {
			if (!frontier.isMemoized(sn.node())) {
				throw new SanityCheckFailureException(String.format("Subnode %s reachable from %s not memoized. [%s, maxVisitedEndPosition=%d, nextPosition=%d]",
						sn, node, frontier, maxVisitedEndPosition, nextPosition()).replace('\n', ' '));
			}
			
			if (sn.node().isReference() && node.node.isReference()) {
				// drop reference-reference transitions as we're only
				// traversing non-reference paths
			} else {
				TraversalNode tn = new TraversalNode(node, sn, sn.isReference() ? ANCHORED_SCORE - sn.weight() : 0);
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
			// only set up starting paths for nodes that should be memoized
			if (frontier.isMemoized(child)) {
				TraversalNode tn = new TraversalNode(new KmerPathSubnode(child), child.isReference() ? ANCHORED_SCORE - child.weight() : 0);
				frontier.memoize(tn);
			}
		}
		if (Defaults.SANITY_CHECK_MEMOIZATION) {
			sanityCheckAreRemovedFromPaths(ImmutableSet.of(node));
			sanityCheck();
		}
	}
	/**
	 * Removes a node from the graph.
	 * 
	 * @param node node to removed
	 */
	@Override
	public void remove(Set<KmerPathNode> nodes) {
		frontier.remove(nodes);
		frontier.tracking_lastRemoval().pathsRestarted = restartChildren(nodes);
		if (Defaults.SANITY_CHECK_MEMOIZATION) {
			sanityCheckAreRemovedFromPaths(nodes);
			sanityCheck();
		}
	}
	private int restartChildren(Set<KmerPathNode> nodes) {
		int count = 0;
		for (KmerPathNode node : nodes) {
			// flag the successors as potential starting nodes
			for (KmerPathNode child : node.next()) {
				// only set up starting paths for nodes that should be memoized
				if (nodes.contains(child)) continue;
				if (!frontier.isMemoized(child)) continue;
				count++;
				TraversalNode tn = new TraversalNode(new KmerPathSubnode(child), child.isReference() ? ANCHORED_SCORE - child.weight() : 0);
				frontier.memoize(tn);
			}
		}
		// technically this is a count of memoization calls made, not
		// of the number of nodes that a restart occurred at
		return count;
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
			int loadBefore = nextPosition() + EVIDENCE_WIDTH_LOADING_MULTIPLE * maxEvidenceWidth;
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
		assert(contig.stream().allMatch(sn -> sn.node().isValid()));
		return contig;
	}
	@Override
	public int tracking_memoizedNodeCount() {
		return frontier.tracking_memoizedNodeCount();
	}
	@Override
	public int tracking_frontierSize() {
		return frontier.tracking_frontierSize();
	}
	@Override
	public MemoizationStats tracking_lastRemoval() {
		return frontier.tracking_lastRemoval();
	}
	@Override
	public void exportState(File file) throws IOException {
		frontier.export(file);
	}
	public boolean sanityCheck(Set<KmerPathNode> loadedGraph) {
		for (KmerPathNode node : loadedGraph) {
			if (!frontier.isMemoized(node)) {
				try {
					File f = File.createTempFile("gridss.sanityfailure.memoization.", ".csv");
					//File f2 = File.createTempFile("gridss.sanityfailure.", ".csv");
					log.error(String.format("Sanity check failure. Unmemoized node %s. Dumping memoization to %s", node, f));
					exportState(f);
				} catch (IOException e) {
				}
			}
			assert(frontier.isMemoized(node));
		}
		return (sanityCheck());
	}
	public boolean sanityCheck() {
		for (TraversalNode tn : contigByScore) {
			assert(frontier.memoized(tn.node.node()).contains(tn));
			if (tn.parent == null && tn.node.isReference()) {
				assert(tn.score == ANCHORED_SCORE);
			}
			if (tn.parent != null) {
				TraversalNode parent = tn.parent;
				assert(parent.node.node().isValid());
				assert(frontier.isMemoized(parent.node.node()));
				if (!parent.node.node().isReference()) {
					// frontier goes to the expected parent node
					assert(frontier.memoized(parent.node.node()).contains(parent));
				}
			}
		}
		for (TraversalNode tn : frontierByPathStart) {
			assert(frontier.memoized(tn.node.node()).contains(tn));
		}
		return frontier.sanityCheck();
	}
	/*
	 * Checks that all possible frontier nodes have indeed been traversed 
	 */
	public boolean sanityCheckFrontier() {
		for (TraversalNode tn : frontierByPathStart) {
			assert(tn.node.lastEnd() + 1 >= nextPosition());
		}
		return true;
	}
	/**
	 * Checks that nothing in the graph contains any of the given nodes
	 * @param nodes
	 */
	private boolean sanityCheckAreRemovedFromPaths(Set<KmerPathNode> nodes) {
		Set<TraversalNode> processed = new HashSet<>();
		for (TraversalNode tn : contigByScore) {
			sanityCheckAreRemovedFromPaths_node(nodes, processed, tn);
		}
		for (TraversalNode tn : frontierByPathStart) {
			sanityCheckAreRemovedFromPaths_node(nodes, processed, tn);
		}
		return true;
	}
	private void sanityCheckAreRemovedFromPaths_node(Set<KmerPathNode> excluded, Set<TraversalNode> processed, TraversalNode node) {
		if (processed.contains(node)) return;
		processed.add(node);
		if (excluded.contains(node.node.node())) {
			try {
				File f = File.createTempFile("gridss.sanityfailure.memoization.removal.", ".csv");
				log.error(String.format("Sanity check failure. %s not removed. Dumping memoization to %s", node, f));
				exportState(f);
			} catch (IOException e) {
			}
		}
		if (node.parent != null) {
			sanityCheckAreRemovedFromPaths_node(excluded, processed, node);
		}
	}
	/**
	 * Checks that the two memoizations match. Nodes are considered
	 * equal if they have the same score over the same interval. Since
	 * node traversal order is nondeterministic, the two memoizations
	 * could have different parent paths.
	 * 
	 * Note that TraversalNode counts for a KmerPathNode do not
	 * necessarily match since {[1,2]w=3} and {[1,1]w=3, [2,2]w=3}
	 * are equivalent.  
	 * 
	 * @param caller caller to compare to
	 */
	public void sanityCheckMatches(MemoizedContigCaller caller) {
		NavigableSet<TraversalNode> set1 = new TreeSet<>(TraversalNode.ByKmerScoreStartEnd);
		NavigableSet<TraversalNode> set2 = new TreeSet<>(TraversalNode.ByKmerScoreStartEnd);
		set1.addAll(contigByScore);
		set2.addAll(caller.contigByScore);
		sanityCheckMatches(set1, set2);
		set1 = new TreeSet<>(TraversalNode.ByKmerScoreStartEnd);
		set2 = new TreeSet<>(TraversalNode.ByKmerScoreStartEnd);
		set1.addAll(frontierByPathStart);
		set2.addAll(caller.frontierByPathStart);
		sanityCheckMatches(set1, set2);
	}
	public static void sanityCheckMatches(NavigableSet<TraversalNode> set1, NavigableSet<TraversalNode> set2) {
		for (TraversalNode tn : set1) {
			sanityCheckContains(tn, set2);
		}
		for (TraversalNode tn : set2) {
			sanityCheckContains(tn, set1);
		}
	}
	public static void sanityCheckContains(TraversalNode node, NavigableSet<TraversalNode> set) {
		Iterator<TraversalNode> it = null;
		TraversalNode start = set.floor(node);
		if (start != null) {
			it = set.tailSet(start).iterator();
		} else {
			it = set.iterator();
		}
		int overlap = 0;
		while (it.hasNext()) {
			TraversalNode tn = it.next();
			if (tn.node.firstKmer() < node.node.firstKmer()) continue;
			if (tn.node.firstKmer() > node.node.firstKmer()) break;
			if (tn.node.firstEnd() < node.node.firstStart()) continue;
			if (tn.node.firstStart() > node.node.firstEnd()) continue;
			if (tn.node.weight() == node.node.weight()) {
				overlap += IntervalUtil.overlapsWidthClosed(tn.node.firstStart(), tn.node.firstEnd(), node.node.firstStart(), node.node.firstEnd()); 
			}
		}
		if (overlap != node.node.firstEnd() - node.node.firstStart() + 1) {
			log.debug(String.format("No matching memoized path of weight %d for %s", node.node.weight(), node.node));
		}
		assert(overlap == node.node.firstEnd() - node.node.firstStart() + 1);
	}
}
