package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeSet;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * Memoizes longest weighted path
 * 
 * Unlike BestNonReferenceContigCaller, memoization persists
 * across contig calls.
 * 
 * TODO: reduce descendant recalculation by using the previously
 * memoized score if all previous nodes are marked as unchanged
 * (either by not being flushed, or by recalculation resulting
 * in no change.
 * 
 * @author Daniel Cameron
 *
 */
public class PathCache {
	/**
	 * Memoized path scores in order of descending score
	 */
	private SortedSet<TraversalNode> byScore = new TreeSet<>(TraversalNode.ByScoreDescPosition);
	/**
	 * Memoized path scores for each node
	 */
	private Multimap<KmerPathNode, TraversalNode> memoize = HashMultimap.create();
	/**
	 * Nodes which have not yet been memoized
	 */
	private Set<KmerPathNode> unmemoized = new HashSet<>();
	private Set<KmerPathNode> descendantsFlushed = new HashSet<>();
	private int earliestIncompleteContigStartPosition;
	private int incompletePathsEndAfterPosition;
	/**
	 * Removes the given node from the graph
	 */
	public void remove(KmerPathNode node) {
		flushNode(node);
		flushParents(node);
		flushDescendants(node);
		unmemoized.remove(node);
	}
	/**
	 * Adds a new node to the graph
	 * @param node
	 */
	public void add(KmerPathNode node) {
		flushNode(node);
		flushParents(node);
		flushDescendants(node);
	}
	/**
	 * Flushes the memoized paths for the given node
	 * @param node node to unmemoize scores for
	 * @return true if memoized paths were removed, false if no paths were memoized
	 */
	private boolean flushNode(KmerPathNode node) {
		if (!node.isReference()) {
			// Remove from memoization and high scores
			return byScore.removeAll(memoize.removeAll(node));
		}
		unmemoized.add(node);
		return false;
	}
	private void flushParents(KmerPathNode node) {
		if (node.isReference()) {
			for (KmerPathNode prev : node.prev()) {
				if (!prev.isReference()) {
					// path leading to node could now be anchored/unanchored
					// so we need to rescore it
					flushNode(prev);
				}
			}
		}
	}
	private void flushDescendants(KmerPathNode node) {
		Stack<KmerPathNode> frontier = new Stack<KmerPathNode>();
		frontier.add(node);
		while (!frontier.isEmpty()) {
			KmerPathNode n = frontier.pop();
			if (!descendantsFlushed.contains(n)) {
				descendantsFlushed.add(n);
				for (KmerPathNode next : n.next()) {
					if (!next.isReference()) {
						flushNode(next);
						frontier.add(next);
					}
				}
			} else {
				assert(!memoize.containsKey(n));
			}
		}
	}
	/**
	 * Returns the best completed contig encountered in the graph
	 * 
	 * Calls optimal breakend contigs from a subgraph
	 * 
	 *  ---- best terminal
	 *     ****** overlapping read (maxEvidenceWidth)
	 *          -------------------- non-terminal path
	 *                              | loadedBeforePosition
	 *          ^- earliest non-terminal path
	 *
	 * Can only call best terminal when
	 * terminal.lastEnd() + maxEvidenceWidth < earliest non-terminal firstStart()
	 * 
	 * @param unprocessedPosition starting position of first node
	 * 	that has not been loaded into the graph  
	 * @param maxEvidenceDistance maximum distance between start position of the first kmer
	 *  of a read and the end position of the last kmer of a read 
	 * @return highest scoring contig, null if no contigs exist, or
	 *  the highest scoring contig could have evidence overlap with
	 *  a contig not yet fully traversed.
	 */
	public Contig bestContig(int unprocessedPosition, int maxEvidenceDistance) {
		earliestIncompleteContigStartPosition = unprocessedPosition;
		incompletePathsEndAfterPosition = unprocessedPosition - maxEvidenceDistance - 1;
		descendantsFlushed.clear();
		recalculate();
		if (byScore.isEmpty()) {
			return null;
		}
		TraversalNode best = byScore.first();
		if (best.node.lastEnd() >= earliestIncompleteContigStartPosition - maxEvidenceDistance) {
			// best contig could overlap a contig that has not yet
			// been fully traversed -> best contig not necessarily
			// the best contig for the supporting evidence
			return null;
		}
		return new Contig(byScore.first());
	}
	/**
	 * Repopulates the memoization lookup 
	 */
	private void recalculate() {
		// load all starting paths
		
		// transform list of nodes to sorted frontier of paths to process
			// all new starting paths
			// all memoized paths in prev nodes
		
		// Advance frontier (reusing BestNonReferenceContigCaller logic)
		PriorityQueue<TraversalNode> frontier; // head of queue = earliest start position
		while (!frontier.isEmpty()) {
			TraversalNode tn = frontier.poll();
			for (TraversalNode n : sucessors(tn)) {
				// memoize as per BestNonReferenceContigCaller
				// if memoization updated
					// TODO: do we need to special case terminal nodes?
					frontier.add(n);
			}
		}
	}
}
