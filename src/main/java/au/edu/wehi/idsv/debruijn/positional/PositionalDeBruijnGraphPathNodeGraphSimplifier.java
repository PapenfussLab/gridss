package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntRBTreeSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.SortedSet;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNodeUtil;
import au.edu.wehi.idsv.debruijn.HighestWeightSimilarPath;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.AssemblyPaperFigure;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
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
 * Minimal output representation:
 * Both input and output graphs has the minimum number of KmerPathNodes required to represent the graph (invariant A)
 * 
 * Without loss of generality, consider merging two nodes N_1, N_2 with N_1 \in next(M): A holds before merge
 *  - merging nodes N_1 and N_2 into N' will cause A to be violated iff:
 *   - adjacent node M can be collapsed into N', but could not be collapsed into N_1.
 * To merge nodes:
 *  start(N_1) = start(N_2)
 *  end(N_1) = end(N_2)
 *  kmers with hamming distance threshold
 *  |next(N_12)| >= |next(N_1)| degree does not decrease
 *  M degree reduced by 1 iff N_2 \in next(M)
 *  	N_12 can be merged with M 
 *  
 *  Reducing with subsequence kmers is (TODO: formal proof) a local operation
 *  but for A to hold, nodes must also be merged with matching kmers in flanking positions:
 *  
 *  Consider the following graph:
 *  [11-15] 1 GTAC (node A)
 *  [16-20] 2 GTAC (node B)
 *  [11-15] 1 GGAC (node C)
 *  [10-19] 2 GGTA (node D)
 *  
 *  Presume graph reduction merges C into A:
 *  [11-15] 2 GTAC (node A,C)
 *  [16-20] 2 GTAC (node B)
 *  [10-19] 2 GGTA (node D)
 *  This violates invariant A.
 *  
 *  Merging AC with B (since they differ only in validity intervals, and these intervals are adjacent) results in:  
 *  [11-20] 2 GTAC (node A,B,C)
 *  [10-19] 2 GGTA (node D)
 *  This still violates invariant A.
 *  
 *  Merging ABC with D:
 *  [10-19] 4 GGTAC (nodes A,B,C,D)
 *  
 * A larger example would result in additional nodes being further merged.
 * The number of possible merges is unbounded, but since each merge results in a node containing N_1
 * 
 * @author Daniel Cameron
 *
 */
public class PositionalDeBruijnGraphPathNodeGraphSimplifier implements Iterator<KmerPathNode>, DeBruijnGraph<KmerPathSubnode> {
	private final PeekingIterator<KmerPathNode> underlying;
	private final AssemblyParameters ap;
	private final int maxSupportWidth;
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
	public PositionalDeBruijnGraphPathNodeGraphSimplifier(
			Iterator<KmerPathNode> it,
			AssemblyParameters ap,
			int maxSupportWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.ap = ap;
		this.maxSupportWidth = maxSupportWidth;
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
	private int unchangedBeforePosition() {
		int earliestCollapseEndPosition = currentPosition - ap.positionalMaxPathCollapseLengthInBases(maxSupportWidth);
		int earliestCollapseStartPosition = earliestCollapseEndPosition - maxSupportWidth;
		// collapse can cause a cascade of nodes to merge together (see class-level docs) 
		int earliestMergeStartPosition = earliestCollapseStartPosition - ap.positionalMaxPathLengthInBases(maxSupportWidth);
		return earliestMergeStartPosition - 1;
	}
	private void ensureBuffer() {
		while (currentPosition < Integer.MAX_VALUE && processed.isEmpty() || processed.peek().endPosition() > unchangedBeforePosition()) {
			// load nodes into graph until we are guaranteed
			// to be able to traverse all branches/leaves
			// rooted at the first node
			if (underlying.hasNext()) {
				currentPosition = underlying.peek().startPosition(0);
				addToGraph();
			} else {
				currentPosition = Integer.MAX_VALUE;
			}
			while (collapse() > 0) {
				// collapse as much as we can
			}
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
			collapseCount += collapse(node, false, false, true);
		}
		return collapseCount;
	}
	private int collapseForwardLeaves() {
		int collapseCount = 0;
		while (!activeLeafForward.isEmpty() && activeLeafForward.peek().forwardLeafProcessingPosition <= currentPosition) {
			UnprocessedNode node = activeLeafForward.poll();
			collapseCount += collapse(node, true, false, false);
		}
		return collapseCount;
	}
	private int collapseBackwardLeaves() {
		int collapseCount = 0;
		while (!activeBranchEnd.isEmpty() && activeBranchEnd.peek().backwardLeafProcessingPosition <= currentPosition) {
			UnprocessedNode node = activeLeafBackward.poll();
			collapseCount += collapse(node, false, true, false);
		}
		return collapseCount;
	}
	private int collapse(UnprocessedNode node, boolean collapseForwardLeaf, boolean collapseBackwardLeaf, boolean collapseBranches) {
		if (!node.node.isValid()) return 0;
		if (!node.isUnchanged()) {
			addToQueues(node.node);
			return 0;
		}
		int collapseCount = 0;
		while (tryCollapse(node, collapseForwardLeaf, collapseBackwardLeaf, collapseBranches)) {
			collapseCount++;
		}
		return collapseCount;
	}
	private boolean tryCollapse(UnprocessedNode node, boolean collapseForwardLeaf, boolean collapseBackwardLeaf, boolean collapseBranches) {
		if (!node.node.isValid()) return false;
		EdgeIntervalIterator eii = new EdgeIntervalIterator(node.node);
		while (eii.isValid()) {
			if (collapseForwardLeaf && eii.currentNext().size() == 0 && eii.currentPrev().size() == 1) {
				assert(node.forwardLeafProcessingPosition <= currentPosition);
				if (tryCollapseForwardLeaf(node.node, eii.currentFirstKmerStartPosition(), eii.currentFirstKmerEndPosition())) {
					return true;
				}
			} else if (collapseBackwardLeaf && eii.currentNext().size() == 1 && eii.currentPrev().size() == 0) {
				assert(node.backwardLeafProcessingPosition <= currentPosition);
				if (tryCollapseBackwardLeaf(node.node, eii.currentNext().iterator().next(), eii.currentFirstKmerStartPosition(), eii.currentFirstKmerEndPosition())) {
					return true;
				}
			} else if (eii.currentPrev().size() > 1) {
				assert(node.branchEndProcessingPosition <= currentPosition);
				if (tryCollapseBranch(node.node, eii.currentFirstKmerStartPosition(), eii.currentFirstKmerEndPosition())) {
					return true;
				}
			}
			eii.advance();
		}
		return false;
	}
	private boolean tryCollapseForwardLeaf(KmerPathNode node, int firstKmerStart, int firstKmerEnd) {
		KmerPathSubnode leaf = new KmerPathSubnode(node, firstKmerStart, firstKmerEnd);
		List<KmerPathSubnode> collapsePath = new HighestWeightSimilarPath<KmerPathSubnode>(ap.maxBaseMismatchForCollapse, leaf, true, this, null).find();
		if (collapsePath != null) {
			merge(ImmutableList.of(leaf), collapsePath, 0, 0);
			return true;
		}
		return false;
	}
	/**
	 * Merges the given source path into the target path 
	 * @param sourcePath path to merge
	 * @param targetPath merge destination
	 * @param sourceSkipKmers number of starting kmers to ignore in source
	 * @param targetSkipKmers number of starting kmers to ignore in target
	 */
	private void merge(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath, int sourceSkipKmers, int targetSkipKmers) {
		// TODO
		throw new RuntimeException("NYI");
	}
	private void merge(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath) {
		List<KmerPathNode> source = positionSplit(sourcePath);
		List<KmerPathNode> target = positionSplit(targetPath);
		IntSortedSet kmerStartPositions = new IntRBTreeSet();
		for (KmerPathNode n : source) {
			kmerStartPositions.add(n.startPosition(0));
		}
		for (KmerPathNode n : target) {
			kmerStartPositions.add(n.startPosition(0));
		}
		source = lengthSplit(kmerStartPositions, source);
		target = lengthSplit(kmerStartPositions, target);
		
		// merge the remainders together
		assert(source.size() == target.size());
		for (int i = 0; i < source.size(); i++) {
			KmerPathNode toMerge = source.get(i);
			KmerPathNode into = target.get(i);
			merge(toMerge, into);
			
		}
	}
	private List<KmerPathNode> lengthSplit(IntSortedSet startPositions, List<KmerPathNode> path) {
		List<KmerPathNode> result = new ArrayList<KmerPathNode>(startPositions.size());
		Iterator<KmerPathNode> it = path.iterator();
		assert(it.hasNext());
		KmerPathNode current = it.next();
		for (int startPosition : startPositions) {
			// if this happens, we overran our offset, or the path starts at the
			// wrong place. Both should never happen
			assert(current.startPosition(0) <= startPosition);
			while (current.endPosition(0) < startPosition) {
				result.add(current);
				assert(it.hasNext());
				current = it.next();
			}
			if (current.startPosition(0) == startPosition) {
				// nothing to do, we already start at the given position
				continue;
			}
			// need to split
			assert(current.endPosition(0) <= startPosition);
			KmerPathNode split = current.splitAtLength(startPosition - current.startPosition(0));
			result.add(split);
			assert(current.startPosition(0) == startPosition);
		}
		result.add(current);
		while (it.hasNext()) {
			result.add(it.next());
		}
		return result;
	}
	private List<KmerPathNode> positionSplit(List<KmerPathSubnode> path) {
		List<KmerPathNode> list = new ArrayList<KmerPathNode>(path.size());
		for (KmerPathSubnode n : path) {
			list.add(positionSplit(n));
		}
		return list;
	}
	/**
	 * Splits the containing KmerPathNode along the given lines
	 * @param n
	 * @return
	 */
	private KmerPathNode positionSplit(KmerPathSubnode n) {
		KmerPathNode pn = n.node();
		int preLength = n.firstKmerStartPosition() - pn.startPosition(0);
		if (preLength != 0) {
			KmerPathNode preNode = pn.splitAtLength(preLength);
			addToQueues(preNode);
		}
		int postLength = pn.endPosition(0) - n.firstKmerEndPosition();
		if (postLength != 0) {
			pn = pn.splitAtLength(preLength);
			addToQueues(pn);
		}
		return pn;
	}
	private void merge(KmerPathNode toMerge, KmerPathNode into) {
		assert(toMerge.startPosition() == into.startPosition());
		assert(toMerge.endPosition() == into.endPosition());
		assert(toMerge.length() == into.length());
		// merge nodes
		// TODO into.include(toMerge);
		
		// merge neighbours
		// TODO shrinkNeighbours();
	}
	private boolean tryCollapseBackwardLeaf(KmerPathNode node, KmerPathNode anchor, int firstKmerStart, int firstKmerEnd) {
		throw new RuntimeException("NYI");
	}
	private boolean tryCollapseBranch(KmerPathNode node, int firstKmerStart, int firstKmerEnd) {
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
		/**
		 * Node has remained unchanged created
		 * @return
		 */
		public boolean isUnchanged() {
			return version == node.version();
		}
		public final KmerPathNode node;
		public int branchEndProcessingPosition;
		public int forwardLeafProcessingPosition;
		public int backwardLeafProcessingPosition;
		public int processingCompletePosition;
		public int version;
		public UnprocessedNode(KmerPathNode node) {
			this.node = node;
			this.branchEndProcessingPosition = calculateBranchEndProcessingPosition();
			this.forwardLeafProcessingPosition = calculateForwardLeafProcessingPosition();
			this.backwardLeafProcessingPosition = calculateBackwardLeafProcessingPosition();
			this.processingCompletePosition = calcualteProcessingCompletePosition();
			this.version = node.version();
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
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	@Override
	public int getWeight(KmerPathSubnode node) {
		return node.weight();
	}
	@Override
	public List<KmerPathSubnode> next(KmerPathSubnode node) {
		return node.next();
	}
	@Override
	public List<KmerPathSubnode> prev(KmerPathSubnode node) {
		return node.prev();
	}
	@Override
	public void removeNode(KmerPathSubnode node) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void removeEdge(KmerPathSubnode source, KmerPathSubnode sink) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void addNode(KmerPathSubnode node) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void addEdge(KmerPathSubnode source, KmerPathSubnode sink) {
		throw new UnsupportedOperationException();
	}
	@Override
	public Collection<KmerPathSubnode> allNodes() {
		throw new UnsupportedOperationException();
	}
	@Override
	public String toString(Iterable<? extends KmerPathSubnode> path) {
		return String.format("[%d-%d] %s",
			path.iterator().next().firstKmerStartPosition(),
			path.iterator().next().firstKmerEndPosition(),
			KmerEncodingHelper.baseCalls(DeBruijnSequenceGraphNodeUtil.asKmers(path), ap.k));
	}
	@Override
	public int getK() {
		return ap.k;
	}
	@Override
	public long getKmer(KmerPathSubnode node) {
		return node.node().kmer();
	}
	@Override
	public boolean isReference(KmerPathSubnode node) {
		return node.node().isReference();
	}
}