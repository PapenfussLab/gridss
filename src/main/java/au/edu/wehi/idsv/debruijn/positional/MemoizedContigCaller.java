package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.Int2ObjectRBTreeMap;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.SortedSet;
import java.util.TreeSet;

import com.google.common.collect.RangeSet;


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
	 * Memoized path scores in order of descending score
	 */
	private final SortedSet<TraversalNode> byScore = new TreeSet<>(TraversalNode.ByScoreDescPosition);
	private final MemoizedTraverse frontier = new MemoizedTraverse();
	private final Map<KmerPathNode, RangeSet<Integer>> recalculateQueue = new HashMap<>();
	private Collection<TraversalNode> getBestPathsAtNode(KmerPathNode node) { return null; }
	private void queueRecalculate(KmerPathNode node, int firstStart, int firstEnd) { }
	/**
	 * Removes the given traversal node from the memoization
	 * @param tn path to remove
	 * @return true if the node was a highest scoring node, false otherwise
	 */
	private boolean remove(TraversalNode tn) {
		byScore.remove(tn);
		memoized
	}
	private void removePathsAtNode(KmerPathNode node) {
		// remove all paths involving the given node from the graph  
	}
	
	public MemoizedContigCaller(
			Iterator<KmerPathNode> it,
			int maxEvidenceWidth) {
		super(it, maxEvidenceWidth);
	}
	/**
	 * Adds a new node to the graph.
	 * 
	 * New paths 
	 * 
	 * Parent nodes predecessors need to be added to frontier
	 * so paths to the added node can be calculated
	 * 
	 * 
	 * @param node node to removed
	 */
	private void add(KmerPathNode node) {
		queueRecalculate(node, node.firstStart(), node.firstEnd());
	}
	/**
	 * Removes a node from the graph.
	 *  
	 *  
	 * Memoized paths involving the removed node must be removed,
	 * possibly resulting in a new starting path node.
	 * When a memoized path is removed, the best path over the
	 * interval in which that path was the best must be recalculated.
	 * This can be done by adding all alternate paths overlapping
	 * the removed path to the frontier. 
	 * 
	 * @param node node to removed
	 */
	private void remove(KmerPathNode node) {
		for (TraversalNode tn : getBestPathsAtNode(node)) {
			remove(tn);
		}
		// remove descendant node containing this as ancestor
		for (KmerPathNode next : node.next()) {
			unmemoize(next, node);
		}
		// Delete node from all data structures including:
	}
	/**
	 * Unmemoizes all paths through the given node from the given parent
	 * @param node
	 * @param parent
	 */
	private void unmemoize(KmerPathNode node, KmerPathNode parent) {
		for (TraversalNode tn : getBestPathsAtNode(node)) {
			if (tn.parent.node.node() == parent) {
				boolean wasBest = remove(tn);
				if (wasBest) {
					queueRecalculate(node, tn.node.firstStart(), tn.node.firstEnd());
					unmemoize(node, node);
				}
			}
		}
	}
	/**
	 * Attempts to advance the frontier.
	 * 
	 */
	private void advance() {
		PriorityQueue<TraversalNode> frontier = new PriorityQueue<>(KmerNodeUtil.ByWeightDesc);
	}
	@Override
	public ArrayDeque<KmerPathSubnode> bestContig() {
		return null;
	}
}
