package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.SanityCheckFailureException;
import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.TraversalNodeByPathFirstStartEndSubnodeSortedSet;
import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.TraversalNodeByScoreDescPathFirstIdentity;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker.MemoizationStats;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.util.Log;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;


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
	/**
	 * Path scores in order of descending score
	 */
	private final SortedSet<TraversalNode> contigByScore = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new TraversalNodeByScoreDescPathFirstIdentity() : new TreeSet<>(TraversalNode.ByScoreDescPathFirstEndSubnode);
	// We could convert this into an Int2IntSortedMap if we changed MemoizedContigTraverse
	// to only call onFrontierRemove() on nodes that are actually in the frontier
	private final SortedSet<TraversalNode> frontierByPathStart = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new TraversalNodeByPathFirstStartEndSubnodeSortedSet(16) : new TreeSet<>(TraversalNode.ByPathFirstStartScoreEndSubnode);
	private final MemoizedContigTraverse frontier = new MemoizedContigTraverse();
	
	private int contigByScoreBeforePosition_startPosition = Integer.MIN_VALUE;
	private SortedSet<TraversalNode> contigByScoreBeforePosition = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new TraversalNodeByScoreDescPathFirstIdentity() : new TreeSet<>(TraversalNode.ByScoreDescPathFirstEndSubnode);
	/**
	 * Scoring bonus for anchoring the start/end of a contig at a reference node. 
	 */
	private final int anchoredScore;
	private int maxVisitedEndPosition = Integer.MIN_VALUE;
	private class MemoizedContigTraverse extends MemoizedTraverse {
		@Override
		protected void onMemoizeAdd(TraversalNode tn) {
			if (tn.node.isReference()) {
				if (tn.parent == null) {
					// don't track reference-only paths
					return;
				}
				// shouldn't be traversing reference-reference paths
				// as they'll overwrite non-ref > ref paths terminating
				// at the ending reference node
				assert(!tn.parent.node.isReference());
			}
			contigByScore.add(tn);
			if (tn.pathFirstStart() < contigByScoreBeforePosition_startPosition) {
				contigByScoreBeforePosition.add(tn);
			}
		}
		@Override
		protected void onMemoizeRemove(TraversalNode tn) {
			contigByScore.remove(tn);
			if (tn.pathFirstStart() < contigByScoreBeforePosition_startPosition) {
				contigByScoreBeforePosition.remove(tn);
			}
		}
		@Override
		protected void onMemoizeRemove(Collection<TraversalNode> tns) {
			contigByScore.removeAll(tns);
			for (TraversalNode tn : tns) {
				if (tn.pathFirstStart() < contigByScoreBeforePosition_startPosition) {
					contigByScoreBeforePosition.remove(tn);
				}
			}
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
	public MemoizedContigCaller(int anchoredScore, int maxEvidenceSupportIntervalWidth) {
		super(maxEvidenceSupportIntervalWidth);
		this.anchoredScore = anchoredScore;
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
		TraversalNode tn = new TraversalNode(new KmerPathSubnode(node), node.isReference() ? anchoredScore - node.weight() : 0);
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
		if (Defaults.SANITY_CHECK_MEMOIZATION && Defaults.SANITY_CHECK_MEMOIZATION_ALL_OPERATIONS) {
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
	private void advanceFrontier(int unprocessedPosition) {
		// Can only advance frontier if all possible successors are guaranteed to be loaded
		while (!frontier.isEmptyFrontier() && frontier.peekFrontier().node.lastEnd() < unprocessedPosition - 1) {
			TraversalNode tn = frontier.pollFrontier();
			visit(tn, unprocessedPosition);
		}
		if (Defaults.SANITY_CHECK_MEMOIZATION && Defaults.SANITY_CHECK_MEMOIZATION_ALL_OPERATIONS) {
			sanityCheck();
			sanityCheckFrontier(unprocessedPosition);
		}
	}
	/**
	 * Visits the given frontier node.
	 *
	 * @param toVisit node to visit
	 */
	private void visit(TraversalNode toVisit, int unprocessedPosition) {
		TraversalNode node = toVisit;
		if (node.node.isReference()) {
			// Reset reference traversal to a starting path
			node = new TraversalNode(node.node, anchoredScore - node.node.node().weight());
		}
		assert(node.node.lastEnd() + 1 < unprocessedPosition); // successors must be fully defined
		for (KmerPathSubnode sn : node.node.next()) {
			if (!frontier.isMemoized(sn.node())) {
				throw new SanityCheckFailureException(String.format("Subnode %s reachable from %s not memoized. [%s, maxVisitedEndPosition=%d, nextPosition=%d]",
						sn, node, frontier, maxVisitedEndPosition, unprocessedPosition).replace('\n', ' '));
			}
			
			if (sn.node().isReference() && node.node.isReference()) {
				// drop reference-reference transitions as we're only
				// traversing non-reference paths
			} else {
				TraversalNode tn = new TraversalNode(node, sn, sn.isReference() ? anchoredScore - sn.weight() : 0);
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
	 * @param node node to remove
	 */
	@Override
	public void remove(KmerPathNode node) {
		frontier.remove(node);
		// flag the successors as potential starting nodes
		for (KmerPathNode child : node.next()) {
			// only set up starting paths for nodes that should be memoized
			if (frontier.isMemoized(child)) {
				TraversalNode tn = new TraversalNode(new KmerPathSubnode(child), child.isReference() ? anchoredScore - child.weight() : 0);
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
	 * @param nodes nodes to remove
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
				TraversalNode tn = new TraversalNode(new KmerPathSubnode(child), child.isReference() ? anchoredScore - child.weight() : 0);
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
	 * @param unprocessedPosition 
	 */
	private boolean canCallBestContig(int unprocessedPosition) {
		if (contigByScore.isEmpty()) return false;
		if (!frontierByPathStart.isEmpty()) {
			int frontierPathFirstStart = frontierByPathStart.first().pathFirstStart();
			unprocessedPosition = Math.min(unprocessedPosition, frontierPathFirstStart);
		}
		int bestContigLastEnd = contigByScore.first().node.lastEnd(); 
		return bestContigLastEnd < unprocessedPosition - maxEvidenceSupportIntervalWidth - 1;
	}
	private TraversalNode bestTraversal(int unprocessedPosition) {
		advanceFrontier(unprocessedPosition);
		if (!canCallBestContig(unprocessedPosition)) {
			return null;
		}
		if (contigByScore.isEmpty()) return null;
		TraversalNode tn = contigByScore.first();
		return tn;
	}
	private ArrayDeque<KmerPathSubnode> asUnanchoredPath(TraversalNode tn) {
		if (tn == null) return null;
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
	public ArrayDeque<KmerPathSubnode> bestContig(int unprocessedPosition) {
		TraversalNode tn = bestTraversal(unprocessedPosition);
		if (tn == null) return null;
		return asUnanchoredPath(tn);
	}
	/**
	 * Calls the best contig before the given position
	 * @param unprocessedPosition
	 * @param contigStartsBefore position contig must start before
	 * @return
	 */
	public ArrayDeque<KmerPathSubnode> callBestContigStartingBefore(int unprocessedPosition, int contigStartsBefore) {
		advanceFrontier(unprocessedPosition);
		ensureContigByScoreBeforePosition(contigStartsBefore);
		if (contigByScoreBeforePosition.isEmpty()) return null;
		return asUnanchoredPath(contigByScoreBeforePosition.first());
	}
	private void ensureContigByScoreBeforePosition(int contigStartsBefore) {
		if (contigByScoreBeforePosition_startPosition != contigStartsBefore) {
			contigByScoreBeforePosition = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new TraversalNodeByScoreDescPathFirstIdentity() : new TreeSet<>(TraversalNode.ByScoreDescPathFirstEndSubnode);
			contigByScore.stream()
				.filter(n -> n.pathFirstStart() < contigStartsBefore)
				.collect(Collectors.toCollection(() -> contigByScoreBeforePosition));
			contigByScoreBeforePosition_startPosition = contigStartsBefore;
		}
	}
	/**
	 * Returns the earliest path start still in the frontier
	 * @return
	 */
	public int frontierStart(int unprocessedPosition) {
		advanceFrontier(unprocessedPosition);
		if (frontierByPathStart.isEmpty()) return unprocessedPosition;
		return frontierByPathStart.first().pathFirstStart();
	}
	/**
	 * Returns the longest path still in the frontier 
	 * @param unprocessedPosition
	 * @param startingBefore position frontier path must start before
	 * @return
	 */
	public ArrayDeque<KmerPathSubnode> frontierPath(int unprocessedPosition, int startingBefore) {
		if (!frontierByPathStart.isEmpty() && frontierByPathStart.first().pathFirstStart() < startingBefore) {
			// We could have an early frontier path because we just haven't performed the memoization yet
			advanceFrontier(unprocessedPosition);
			if (!frontierByPathStart.isEmpty() && frontierByPathStart.first().pathFirstStart() < startingBefore) {
				return asUnanchoredPath(frontierByPathStart.first());
			}
		}
		return null;
	}
	@Override
	public int memoizedNodeCount() {
		return frontier.memoizedNodeCount();
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
	public void exportScores(File file) throws IOException {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(file))) {
			writer.write("start,score\n");
			for (TraversalNode tn : contigByScore) {
				writer.write(String.format("%d,%d\n", tn.score, tn.pathFirstStart()));
			}
		}
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
				assert(tn.score == anchoredScore);
			}
			if (tn.parent != null) {
				TraversalNode parent = tn.parent;
				assert(parent.node.node().isValid());
				assert(frontier.isMemoized(parent.node.node()));
				if (!parent.node.node().isReference()) {
					// frontier goes to the expected parent node
					assert(frontier.getParent(tn.node.node(), tn.node.firstStart(), tn.node.firstEnd()) == parent.node.node());
				}
			}
		}
		assert(frontier.tracking_frontierSize() == frontierByPathStart.size());
		for (TraversalNode tn : frontierByPathStart) {
			assert(frontier.memoized(tn.node.node()).contains(tn));
		}
		return frontier.sanityCheck();
	}
	/*
	 * Checks that all possible frontier nodes have indeed been traversed 
	 */
	public boolean sanityCheckFrontier(int unprocessedPosition) {
		for (TraversalNode tn : frontierByPathStart) {
			assert(tn.node.lastEnd() + 1 >= unprocessedPosition);
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
