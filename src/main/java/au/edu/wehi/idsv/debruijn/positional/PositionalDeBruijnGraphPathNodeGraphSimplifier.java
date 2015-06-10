package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;
import com.google.common.primitives.Ints;

/**
 * Graph simplifier that merging similar paths
 * 
 * Input: KmerPathNodes in ascending order of start position of the first kmer
 * 
 * Output: KmerPathNodes in ascending order of start position of the first kmer
 * after graph simplification. Each output node is guaranteed to not be modified
 * by further graph reduction although adjacent node may be modified causing a
 * change in the node edge list after emission.
 *  
 * All input nodes are guaranteed to have all edges defined,
 * but this is not transitive: nodes that have not been returned by the input
 * are not guaranteed to either a) have all edges defined, or b) be of full length
 *  
 * When collapsing 2 path together, all nodes along both paths must be fully defined
 * so the edge list of adjacent nodes can be updated if the node is split.
 * 
 * Let G(pos) be the graph of all nodes such that first_kmer_start_position(n) <= pos
 * 
 * For branch collapse:
 * ... -* - * <- root node we traverse backwards looking for similar paths
 *        /
 * ...   *
 * 
 * Let n be a node such that next(n) = root
 * To collapse a path onto n, n needs to be split along:
 * a) length: alternate path may have nodes of shorter length,
 *     this requires this node to be split into shorter nodes
 * b) start/end: the alternate path may have a smaller positional window of validity,
 *     requiring splitting into multiple validity intervals
 *
 * If we collapse a branch as soon as pos >= first_kmer_end_position of root then:
 * - n in G(pos) since last_kmer_start_position(n) < last_kmer_end_position(n) < first_kmer_start_position(n) <= pos
 * - break n length-wise
 *   - need to prove that the partially constructed nodes adjacent to n will reference the correct split
 *    - yes: neighbours of split before n' are fully defined 
 *    - yes: neighbours of n' are fully defined
 *    - yes: neighbours of post n' are partially defined iff they only connect to post n' -> connect to post n' post-split
 *  
 * Leaf collapse:
 *    * - - - - -   <- leaf     
 *                \    
 * * - * - * - * - * - ?
 *                 ^
 *                root
 * 
 * Conditions for backward leaf collapse match that of path collapse when the
 * root node is considered to be the node adjacent to the leaf
 * 
 *     - - - - - *      <- collapse this leave into a main path
 *   /
 * * - * - * - * - * - ?
 * 
 * Similarly, forward leaf collapse requires fully defined path nodes up to last_kmer_end_position(leaf)
 * 
 * 
 * @author Daniel Cameron
 *
 */
public class PositionalDeBruijnGraphPathNodeGraphSimplifier implements Iterator<KmerPathNode> {
	private final PeekingIterator<KmerPathNode> underlying;
	private final int k;
	private final int maxSupportWidth;
	private final int maxCollapseLength;
	private final int maxBaseMismatchForCollapse;
	private final boolean collapseBubblesOnly;
	private final PriorityQueue<KmerPathNode> processed = new PriorityQueue<KmerPathNode>(1024, KmerPathNode.ByFirstKmerStartPosition);
	private final PriorityQueue<UnprocessedNode> active = new PriorityQueue<UnprocessedNode>(1024, ByProcessingCompletePosition);
	/**
	 * Unprocessed forward branching leaves
	 */
	private final PriorityQueue<UnprocessedNode> activeLeafForward = new PriorityQueue<UnprocessedNode>(1024, ByForwardLeafProcessingPosition);
	/**
	 * Unprocessed end of possible branches.
	 * Since nodes are added in 
	 */
	private final PriorityQueue<UnprocessedNode> activeBranchEnd = new PriorityQueue<UnprocessedNode>(1024, ByBranchEndProcessingPosition);
	/**
	 * Unprocessed backward branching leaves
	 */
	private final PriorityQueue<UnprocessedNode> activeLeafBackward = new PriorityQueue<UnprocessedNode>(1024, ByBackwardLeafProcessingPosition);
	/**
	 * Position at which all nodes starting at or this position have been loaded
	 */
	private int currentPosition = Integer.MIN_VALUE;
	public PositionalDeBruijnGraphPathNodeGraphSimplifier(Iterator<KmerPathNode> it, int k, int maxSupportWidth, int maxCollapseLength, boolean collapseBubblesOnly, int maxBaseMismatchForCollapse) {
		this.underlying = Iterators.peekingIterator(it);
		this.k = k;
		this.maxSupportWidth = maxSupportWidth;
		this.maxCollapseLength = maxCollapseLength;
		this.maxBaseMismatchForCollapse = maxBaseMismatchForCollapse;
		this.collapseBubblesOnly = collapseBubblesOnly;
	}
	@Override
	public boolean hasNext() {
		if (!processed.isEmpty()) return true;
		ensureBuffer();
		return !processed.isEmpty();
	}
	@Override
	public KmerPathNode next() {
		ensureBuffer();
		return processed.poll();
	}
	private void ensureBuffer() {
		while (currentPosition < Integer.MAX_VALUE && processed.isEmpty() || processed.peek().endPosition() > currentPosition - maxCollapseLength - 1
				/* maxSupportWidth required since backward leaves are delayed until the next node is emitted which could be as late as maxSupportWidth
				 * after leaf node is emitted.
				 */
				- maxSupportWidth) {
			// load nodes into graph until we are guaranteed
			// to be able to traverse all branches/leaves
			// rooted at the first node
			if (underlying.hasNext()) {
				currentPosition = underlying.peek().startPosition(0);
				addToGraph();
			} else {
				currentPosition = Integer.MAX_VALUE;
			}
			collapse();
			while (!active.isEmpty() && active.peek().processingCompletePosition <= currentPosition) {
				processed.add(active.poll().node);
			}
		}
	}
	private void addToGraph() {
		while (underlying.hasNext() && underlying.peek().startPosition(0) <= currentPosition) {
			addToQueues(underlying.next());
		}
	}
	private void addToQueues(KmerPathNode node) {
		UnprocessedNode un = new UnprocessedNode(node);
		// can't use more restrictive conditions without
		// fully decoding edge validity intervals
		// with (eg EdgeIntervalIterator)
		// additionally, if nodes are subsequently
		// merge into this, edge counts in both
		// directions can increase
		if (node.prev().size() > 0) {
			// forward leaf = 1 prev, 0 next
			activeLeafForward.add(un);
			// branch end candidate = 2+ prev
			activeBranchEnd.add(un);
		}
		if (node.next().size() > 0) {
			// backward leaf = 0 prev, 1 next
			activeLeafBackward.add(un);
		}
		active.add(un);
	}
	private int collapse() {
		int collapseCount = collapseBranches();
		collapseCount += collapseForwardLeaves();
		collapseCount += collapseBackwardLeaves();
		return collapseCount;
	}
	private int collapseBranches() {
		int collapseCount = 0;
		while (!activeBranchEnd.isEmpty() && activeBranchEnd.peek().branchEndProcessingPosition <= currentPosition) {
			UnprocessedNode node = activeBranchEnd.poll();
			if (tryCollapseBranch(node.node)) {
				collapseCount++;
			}
		}
		return collapseCount;
	}
	private int collapseForwardLeaves() {
		int collapseCount = 0;
		while (!activeLeafForward.isEmpty() && activeLeafForward.peek().forwardLeafProcessingPosition <= currentPosition) {
			UnprocessedNode node = activeLeafForward.poll();
			if (tryCollapseForwardLeaf(node.node)) {
				collapseCount++;
			}
		}
		return collapseCount;
	}
	private int collapseBackwardLeaves() {
		int collapseCount = 0;
		while (!activeBranchEnd.isEmpty() && activeBranchEnd.peek().backwardLeafProcessingPosition <= currentPosition) {
			UnprocessedNode node = activeLeafBackward.poll();
			if (tryCollapseBackwardLeaf(node.node)) {
				collapseCount++;
			}
		}
		return collapseCount;
	}
	private boolean tryCollapseForwardLeaf(KmerPathNode node) {
		if (!node.isValid()) return false;
		// determine the intervals in which this node is a forward leaf
		EdgeIntervalIterator eii = new EdgeIntervalIterator(node);
		IntArrayList starts = new IntArrayList();
		IntArrayList ends = new IntArrayList();
		while (eii.isValid()) {
			if (eii.currentNext().size() == 0 && eii.currentPrev().size() == 1) {
				starts.add(eii.currentFirstKmerStartPosition());
				ends.add(eii.currentFirstKmerEndPosition());
			}
			eii.advance();
		}
		//tryCollapseForwardLeaf(node, start, end);
		throw new RuntimeException("NYI");
	}
	private boolean tryCollapseBackwardLeaf(KmerPathNode node) {
		if (!node.isValid()) return false;
		throw new RuntimeException("NYI");
	}
	private boolean tryCollapseBranch(KmerPathNode node) {
		if (!node.isValid()) return false;
		throw new RuntimeException("NYI");
	}
	/**
	 * Processing positions are the earliest input position
	 * at which the node can be processed.
	 * 
	 * These inital snapshot positions ensure that 
	 * invariants are not void
	 * 
	 * Snapshot values are required to ensure that PriorityQueue invariants
	 * 
	 * @author Daniel Cameron
	 *
	 */
	private class UnprocessedNode {
		public int calculateBranchEndProcessingPosition() {
			return node.endPosition();
		}
		public int calculateForwardLeafProcessingPosition() {
			return node.endPosition(0);
		}
		public int calculateBackwardLeafProcessingPosition() {
			return maxFirstKmerEndPositionOfSuccessor();
		}
		/**
		 * Position at which no more processing will be performed
		 * @return
		 */
		public int calcualteProcessingCompletePosition() {
			return maxFirstKmerEndPositionOfSuccessor();
		}
		private int maxFirstKmerEndPositionOfSuccessor() {
			int pos = Integer.MIN_VALUE;
			for (KmerPathNode pn : node.next()) {
				pos = Math.max(pos, pn.endPosition(0));
			}
			return pos;
		}
		public final KmerPathNode node;
		public int branchEndProcessingPosition;
		public int forwardLeafProcessingPosition;
		public int backwardLeafProcessingPosition;
		public int processingCompletePosition;
		public UnprocessedNode(KmerPathNode node) {
			this.node = node;
			this.branchEndProcessingPosition = calculateBranchEndProcessingPosition();
			this.forwardLeafProcessingPosition = calculateForwardLeafProcessingPosition();
			this.backwardLeafProcessingPosition = calculateBackwardLeafProcessingPosition();
			this.processingCompletePosition = calcualteProcessingCompletePosition();
		}
	}
	private static Ordering<UnprocessedNode> ByBranchEndProcessingPosition = new Ordering<UnprocessedNode>() {
		@Override
		public int compare(UnprocessedNode left, UnprocessedNode right) {
			return Ints.compare(left.branchEndProcessingPosition, right.branchEndProcessingPosition);
		}
	};
	private static Ordering<UnprocessedNode> ByForwardLeafProcessingPosition = new Ordering<UnprocessedNode>() {
		@Override
		public int compare(UnprocessedNode left, UnprocessedNode right) {
			return Ints.compare(left.forwardLeafProcessingPosition, right.forwardLeafProcessingPosition);
		}
	};
	private static Ordering<UnprocessedNode> ByBackwardLeafProcessingPosition = new Ordering<UnprocessedNode>() {
		@Override
		public int compare(UnprocessedNode left, UnprocessedNode right) {
			return Ints.compare(left.backwardLeafProcessingPosition, right.backwardLeafProcessingPosition);
		}
	};
	private static Ordering<UnprocessedNode> ByProcessingCompletePosition = new Ordering<UnprocessedNode>() {
		@Override
		public int compare(UnprocessedNode left, UnprocessedNode right) {
			return Ints.compare(left.processingCompletePosition, right.processingCompletePosition);
		}
	};
}