package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.SortedSet;
import java.util.TreeSet;


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
	/**
	 * Paths that require further calculation. Initially
	 * 
	 * 
	 */
	private final PriorityQueue<KmerPathSubnode> frontier = new PriorityQueue<>(KmerNodeUtil.ByWeightDesc);
	public MemoizedContigCaller(
			Iterator<KmerPathNode> it,
			int maxEvidenceWidth) {
		super(it, maxEvidenceWidth);
	}
	/**
	 * Adds a new node to the graph.
	 * 
	 * The addition of a new node partially
	 * invalidates the memoization.
	 * 
	 * - 
	 * - the best paths from the added node must
	 * be traversed across all descendant nodes 
	 *
	 * 
	 * If the added node is terminal
	 * 
	 * 
	 * @param node
	 */
	private void add(KmerPathNode node) {
		
	}
	@Override
	public ArrayDeque<KmerPathSubnode> bestContig() {
		return null;
	}
}
