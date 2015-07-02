package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntRBTreeSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.Set;
import java.util.TreeSet;

import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNodeUtil;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNodeBasePath.TraversalNode;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

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
 * For simplicity, each node processed for all collapse types once:
 *                             |----maxSupportWidth+maxPathLength----|
 *    |----maxcollapseLength---|                                     |----maxcollapseLength---|
 *    ^                                                              ^                        ^
 *    |                                      ^^^^                 process                     |
 * unchanged                        node being processed          offset                  input position 
 *  offset
 * 
 * Note: the output graph is not minimal and may contain adjacent nodes that could be merged
 * 
 * @author Daniel Cameron
 *
 */
public class PathCollapseIterator implements Iterator<KmerPathNode>, DeBruijnGraph<KmerPathSubnode> {
	private final PeekingIterator<KmerPathNode> underlying;
	private final int k;
	private final int processOffset;
	private final int maxCollapseLength;
	private final int maxBasesMismatch;
	private final boolean bubblesAndLeavesOnly;
	private final NavigableSet<KmerPathNode> processed = new TreeSet<KmerPathNode>(KmerNodeUtil.ByFirstStartEndKmerReference);
	private final NavigableSet<KmerPathNode> unprocessed = new TreeSet<KmerPathNode>(KmerNodeUtil.ByLastEndStartKmerReference );
	private int inputPosition = Integer.MIN_VALUE;
	private int maxNodeWidth = 0;
	private int maxNodeLength = 0;
	private int emitOffset() {
		// records ending before this position cannot be changed by subsequent operations
		int unchangedOffset = processOffset + maxNodeLength + maxNodeWidth + maxCollapseLength + 1;
		// records are added to the emit queue when their last kmer end position position is before unchangedOffset
		// we need to resort so they are output in the correct start kmer order
		return unchangedOffset + maxNodeLength + maxNodeWidth + 1; 
	}
	public PathCollapseIterator(
			Iterator<KmerPathNode> it,
			int k,
			int maxPathCollapseLength,
			int maxBasesMismatch,
			boolean bubblesAndLeavesOnly) {
		this.underlying = Iterators.peekingIterator(it);
		this.k = k;
		this.maxBasesMismatch = maxBasesMismatch;
		this.maxCollapseLength = maxPathCollapseLength;
		this.processOffset = maxCollapseLength + 1;
		this.bubblesAndLeavesOnly = bubblesAndLeavesOnly;
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		boolean hasNext = !processed.isEmpty();
		if (!hasNext) {
			assert(unprocessed.isEmpty());
			assert(!underlying.hasNext());
		}
		return hasNext;
	}
	@Override
	public KmerPathNode next() {
		ensureBuffer();
		KmerPathNode node = processed.pollFirst();
		return node;
	}
	private void ensureBuffer() {
		while (inputPosition < Integer.MAX_VALUE && (processed.isEmpty() || processed.first().firstStart() > inputPosition - emitOffset())) {
			// advance graph position
			if (underlying.hasNext()) {
				inputPosition = underlying.peek().firstStart();
			} else {
				inputPosition = Integer.MAX_VALUE;
			}
			loadGraphNodes();
			while (collapse() > 0) { } // collapse as much as we can
		}
	}
	/**
	 * Adds all nodes with first kmer starting at or before the current inputPosition
	 * to processing queue
	 */
	private void loadGraphNodes() {
		while (underlying.hasNext() && underlying.peek().firstStart() <= inputPosition) {
			KmerPathNode node = underlying.next();
			maxNodeWidth = Math.max(maxNodeWidth, node.width());
			maxNodeLength = Math.max(maxNodeLength, node.length());
			unprocessed.add(node);
		}
	}
	private int collapse() {
		int collapseCount = 0;
		while (!unprocessed.isEmpty() && unprocessed.first().lastEnd() < inputPosition - processOffset) {
			if (collapseNext()) {
				collapseCount++;
			}
		}
		return collapseCount;
	}
	private boolean collapseNext() {
		KmerPathNode node = unprocessed.pollFirst();
		processed.add(node);
		return collapse(node);
	}
	private boolean collapse(KmerPathNode node) {
		KmerPathSubnode root = new KmerPathSubnode(node);
		List<KmerPathSubnode> nextNodes = root.next();
		for (int i = 0; i < nextNodes.size(); i++) {
			for (int j = i + 1; j < nextNodes.size(); j++) {
				if (collapseSimilarPath(node, nextNodes.get(i), nextNodes.get(j), true, true, true)) return true;
			}
		}
		List<KmerPathSubnode> prevNodes = root.next();
		for (int i = 0; i < prevNodes.size(); i++) {
			for (int j = i + 1; j < prevNodes.size(); j++) {
				if (collapseSimilarPath(node, prevNodes.get(i), prevNodes.get(j), true, false, false)) return true;
			}
		}
		return false;
	}
	private boolean collapseSimilarPath(KmerPathNode root, KmerPathSubnode startNodeA, KmerPathSubnode startNodeB, boolean findLeaf, boolean findCommonChild, boolean traverseForward) {
		KmerPathNodePath pathA = new KmerPathNodePath(startNodeA, traverseForward, maxCollapseLength);
		KmerPathNodePath pathB = new KmerPathNodePath(startNodeB, traverseForward, maxCollapseLength);
		if (pathA.pathLength() <= maxCollapseLength && pathB.pathLength() <= maxCollapseLength) {
			Set<KmerPathNode> set = Collections.newSetFromMap(new IdentityHashMap<KmerPathNode, Boolean>());
			set.add(root);
			set.add(pathA.headPath());
			set.add(pathB.headPath());
			if (set.size() == 3) {
				return collapseSimilarPath(set, pathA, pathB, findLeaf, findCommonChild, traverseForward);
			}
		}
		return false;
	}
	/**
	 * Recusive node traversal of path trees looking for similar paths
	 * 
	 * Need to simultaneous traverse across both trees comparing all possible path combinations until a match is found. 
	 * 
	 * Invariant: pathA and pathB are unchanged when returning false
	 * 
	 * Although the underlying graph is a DAG, KmerPathSubnode traversal is not. The same
	 * underlying KmerPathNode can be present multiple times in either path. Example scenario:
	 * [1,10] AAAA -> [2,11] AAAA -> [3,12] AAAG = [1,10] AAAAG 
	 * [1,10] AAAA -> [2,11] AAAG -> [3,12] AAGG = [1,10] AAAGG
	 * only 1 kmer difference so should collapse but when merging, we are unable to fragment
	 * [1,11] AAAA in such a way that each KmerPathSubnode contains a single KmerPathNode
	 * per KmerPathSubnode, we have to fragment the KmerPathNode across all boundaries.
	 * For now, we handle this by not collapsing if we encounter a KmerPathNode repeat. 
	 * 
	 */
	private boolean collapseSimilarPath(Set<KmerPathNode> onPath, KmerPathNodePath pathA, KmerPathNodePath pathB, boolean findLeaf, boolean findCommonChild, boolean traverseForward) {
		// paths don't share any common interval - no way for them to be the same length and still overlap
		if (!pathsOverlap(pathA, pathB)) return false;
		// paths have too many bases different
		if (!areSimilarPartialPaths(pathA, pathB, traverseForward)) return false;
		if (tryCollapse(pathA, pathB, findLeaf, findCommonChild, traverseForward)) return true;
		if (pathA.headPath() == pathB.headPath()) return false; // common child found and no leaf -> no need to traverse further 
		int pathAlength = pathA.pathLength();
		int pathBlength = pathB.pathLength();
		if (pathA.pathLength() <= pathB.pathLength()) {
			while (pathA.dfsNextChild()) {
				if (!onPath.contains(pathA.headPath()) || pathA.headPath() == pathB.headPath()) {
					onPath.add(pathA.headPath());
					pathB.dfsResetChildTraversal();
					if (collapseSimilarPath(onPath, pathA, pathB, findLeaf, findCommonChild, traverseForward)) return true;
					onPath.remove(pathA.headPath());
				}
				pathA.dfsPop();
			}
			//pathA.dfsPop(); // done with this node
		} else {
			while (pathB.dfsNextChild()) {
				if (!onPath.contains(pathB.headPath()) || pathA.headPath() == pathB.headPath()) {
					onPath.add(pathB.headPath());
					if (collapseSimilarPath(onPath, pathA, pathB, findLeaf, findCommonChild, traverseForward)) return true;
					onPath.remove(pathB.headPath());
				}
				pathB.dfsPop();
			}
		}
		assert(pathA.pathLength() == pathAlength);
		assert(pathB.pathLength() == pathBlength);
		return false;
	}
	private boolean areSimilarPartialPaths(KmerPathNodePath pathA, KmerPathNodePath pathB, boolean traverseForward) {
		int basesDifference;
		if (traverseForward) {
			basesDifference = DeBruijnSequenceGraphNodeUtil.basesDifferent(k, pathA.currentPath(), pathB.currentPath()); 
		} else {
			basesDifference = DeBruijnSequenceGraphNodeUtil.reverseBasesDifferent(k, pathA.currentPath(), pathB.currentPath());
			//basesDifference = DeBruijnSequenceGraphNodeUtil.basesDifferent(k,
			//		StreamSupport.stream(Spliterators.spliteratorUnknownSize(pathA.currentPath().descendingIterator(), Spliterator.ORDERED), false)
			//			.flatMapToLong(n -> IntStream.rangeClosed(n.length() - 1, 0).mapToLong(i -> n.kmer(-i))),
			//		StreamSupport.stream(Spliterators.spliteratorUnknownSize(pathA.currentPath().descendingIterator(), Spliterator.ORDERED), false)
			//			.flatMapToLong(n -> IntStream.rangeClosed(n.length() - 1, 0).mapToLong(i -> n.kmer(-i))),
			//		false);
		}
		return basesDifference <= maxBasesMismatch; 
	}
	private boolean pathsOverlap(KmerPathNodePath pathA, KmerPathNodePath pathB) {
		boolean overlaps = IntervalUtil.overlapsClosed(
				pathA.headNode().startPositionOfAnchorKmer(), pathA.headNode().endPositionOfAnchorKmer(),
				pathB.headNode().startPositionOfAnchorKmer(), pathB.headNode().endPositionOfAnchorKmer());
		return overlaps;
	}
	/**
	 * 
	 * @param pathA
	 * @param pathB
	 * @param findLeaf
	 * @param findBubble
	 * @param traverseForward
	 * @return true if paths were collapsed, false otherwise
	 */
	private boolean tryCollapse(KmerPathNodePath pathA, KmerPathNodePath pathB, boolean findLeaf, boolean findCommonChild, boolean traverseForward) {
		assert(findLeaf || findCommonChild);
		if (pathA.headPath() == pathB.headPath()) {
			if (findCommonChild && pathA.pathLength() == pathB.pathLength()) {
				// remove common trailing node
				pathA.dfsPop();
				pathB.dfsPop();
				assert(repeatedKmerPathNodeCount(pathA, pathB) == 0);
				List<KmerPathSubnode> lA = new ArrayList<KmerPathSubnode>(pathA.headNode().overlapping(pathB.headNode()).asSubnodes());
				List<KmerPathSubnode> lB = new ArrayList<KmerPathSubnode>(pathB.headNode().overlapping(pathA.headNode()).asSubnodes());
				if (pathA.pathWeight() < pathB.pathWeight()) {
					if (!bubblesAndLeavesOnly || isBubblePath(lA)) {
						merge(lA, lB, 0, 0);
						return true;
					}
				} else {
					if (!bubblesAndLeavesOnly || isBubblePath(lB)) {
						merge(lB, lA, 0, 0);
						return true;
					}
				}
			}
		} else {
			if (tryLeafCollapse(pathA, pathB, traverseForward)) return true;
			if (tryLeafCollapse(pathB, pathA, traverseForward)) return true;
		}
		return false;
	}
	/**
	 * A path is considered a bubble if it each node only has a single source and successor 
	 * @param path path excluding starting node of bubble, but including ending node
	 * @return true if the path is a bubble, false otherwise
	 */
	private boolean isBubblePath(List<KmerPathSubnode> path) {
		for (int i = 0; i < path.size() - 1; i++) {
			KmerPathSubnode sn = path.get(i);
			if (sn.next().size() != 1) return false;
			if (sn.prev().size() != 1) return false;
		}
		return true;
	}
	private int repeatedKmerPathNodeCount(KmerPathNodePath... paths) {
		HashSet<KmerPathNode> set = new HashSet<KmerPathNode>();
		int nodeCount = 0;
		for (KmerPathNodePath path : paths) {
			nodeCount += path.currentPath().size();
			set.addAll(path.currentPath());
		}
		// we if shrunk in size then we have a repeat
		return nodeCount - set.size();
	}
	private boolean tryLeafCollapse(KmerPathNodePath leaf, KmerPathNodePath path, boolean traverseForward) {
		// leaf can't be longer than the path
		if (leaf.pathLength() > path.pathLength()) return false;
		if (leaf.pathWeight() > path.pathWeight()) return false;
		TraversalNode firstLeaf = leaf.headNode().overlapping(path.headNode()).firstTerminalLeaf();
		if (firstLeaf == null) return false;
		assert(repeatedKmerPathNodeCount(leaf, path) == 0);
		int leafSkip = 0;
		int pathSkip = traverseForward ? 0 : path.pathLength() - leaf.pathLength();
		merge(new ArrayList<KmerPathSubnode>(firstLeaf.asSubnodes()),
				new ArrayList<KmerPathSubnode>(path.headNode().overlapping(firstLeaf).asSubnodes()),
				leafSkip, pathSkip);
		return true;
	}
	/**
	 * Merges the given source path into the target path 
	 * @param sourcePath path to merge
	 * @param targetPath merge destination
	 * @param sourceSkipKmers number of starting kmers to ignore in source
	 * @param targetSkipKmers number of starting kmers to ignore in target
	 */
	private void merge(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath, int sourceSkipKmers, int targetSkipKmers) {
		trimStartKmers(sourcePath, sourceSkipKmers);
		trimStartKmers(targetPath, targetSkipKmers);
		merge(sourcePath, targetPath);
	}
	private void trimStartKmers(List<KmerPathSubnode> path, int kmerCount) {
		assert(kmerCount >= 0);
		if (kmerCount > 0) {
			lengthSplit(kmerCount, path);
			while (kmerCount > 0) {
				kmerCount -= path.get(0).length();
				path.remove(0);
			}
		}
		assert(kmerCount == 0);
	}
	private void merge(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath) {
		assert(sourcePath.get(0).width() == targetPath.get(0).width());
		assert(sourcePath.get(0).firstStart() == targetPath.get(0).firstStart());
		assert(DeBruijnSequenceGraphNodeUtil.basesDifferent(k, sourcePath, targetPath) <= maxBasesMismatch);
		List<KmerPathNode> source = positionSplit(sourcePath);
		List<KmerPathNode> target = positionSplit(targetPath);
		IntSortedSet kmerStartPositions = new IntRBTreeSet();
		for (KmerPathNode n : source) {
			kmerStartPositions.add(n.firstStart());
			kmerStartPositions.add(n.firstStart() + n.length());
		}
		for (KmerPathNode n : target) {
			kmerStartPositions.add(n.firstStart());
			kmerStartPositions.add(n.firstStart() + n.length());
		}
		source = lengthSplit(kmerStartPositions, source);
		target = lengthSplit(kmerStartPositions, target);
		assert(DeBruijnSequenceGraphNodeUtil.basesDifferent(k, source, target) <= maxBasesMismatch);
		// merge the common nodes
		for (int i = 0; i < Math.min(source.size(), target.size()); i++) {
			KmerPathNode toMerge = source.get(i);
			KmerPathNode into = target.get(i);
			merge(toMerge, into);
		}
	}
	/**
	 * Ensures that the 
	 * @param splitAfter splits after the given number of bases
	 * @param path path to split
	 */
	private void lengthSplit(int splitAfter, List<KmerPathSubnode> path) {
		assert(splitAfter > 0);
		if (splitAfter == 0) return;
		int index = 0;
		int length = 0;
		for (KmerPathSubnode n : path) {
			if (length + n.length() == splitAfter) {
				// already a split at the given position
				return;
			} else if (length + n.length() < splitAfter) {
				// advance to next node
				length += n.length();
				index++;
			} else {
				// split the underlying node
				int splitLength = splitAfter - length;
				KmerPathNode split = lengthSplit(n.node(), splitLength);
				path.add(index, new KmerPathSubnode(split, n.firstStart(), n.firstEnd()));
				path.set(index + 1, new KmerPathSubnode(n.node(), n.firstStart() + splitLength, n.firstEnd() + splitLength));
			}		
		}
	}
	private List<KmerPathNode> lengthSplit(IntSortedSet startPositions, List<KmerPathNode> path) {
		List<KmerPathNode> result = new ArrayList<KmerPathNode>(startPositions.size());
		for (int i = 0; i < path.size(); i++) {
			KmerPathNode pn = path.get(i);
			// break node internally
			for (int breakStartPosition : startPositions.subSet(pn.firstStart() + 1, pn.firstStart() + pn.length())) {
				int breakLength = breakStartPosition - pn.firstStart();
				KmerPathNode split = lengthSplit(pn, breakLength);
				result.add(split);
			}
			result.add(pn);
		}
		return result;
	}
	private KmerPathNode lengthSplit(KmerPathNode node, int length) {
		NavigableSet<KmerPathNode> queue = processed.contains(node) ? processed : unprocessed;
		queue.remove(node);
		KmerPathNode split = node.splitAtLength(length);
		queue.add(split);
		queue.add(node);
		return split;
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
		assert(processed.contains(pn) || unprocessed.contains(pn));
		if (n.firstStart() != pn.firstStart()) {
			NavigableSet<KmerPathNode> queue = processed.contains(pn) ? processed : unprocessed;
			queue.remove(pn);
			KmerPathNode preNode = pn.splitAtStartPosition(n.firstStart());
			queue.add(pn);
			queue.add(preNode);
		}
		if (pn.firstEnd() != n.firstEnd()) {
			NavigableSet<KmerPathNode> queue = processed.contains(pn) ? processed : unprocessed;
			queue.remove(pn);
			KmerPathNode midNode = pn.splitAtStartPosition(n.firstEnd() + 1);
			queue.add(pn);
			queue.add(midNode);
			pn = midNode;
		}
		assert(pn.firstStart() == n.firstStart());
		assert(pn.firstEnd() == n.firstEnd());
		return pn;
	}
	private void merge(KmerPathNode toMerge, KmerPathNode into) {
		assert(toMerge != into);
		assert(toMerge.lastStart() == into.lastStart());
		assert(toMerge.lastEnd() == into.lastEnd());
		assert(toMerge.length() == into.length());
		processed.remove(toMerge);
		unprocessed.remove(toMerge);
		// merge nodes
		into.merge(toMerge);
	}
	public void sanityCheck() {
		assert(false);
		// TODO: (kmer, start) should be unique
		// (kmer, interval) should be non-overlapping
		
		// processed and unprocessed should be mutually exclusive
	}
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
			path.iterator().next().firstStart(),
			path.iterator().next().firstEnd(),
			KmerEncodingHelper.baseCalls(DeBruijnSequenceGraphNodeUtil.asKmers(path), k));
	}
	@Override
	public int getK() {
		return k;
	}
	@Override
	public long getKmer(KmerPathSubnode node) {
		return node.node().lastKmer();
	}
	@Override
	public boolean isReference(KmerPathSubnode node) {
		return node.node().isReference();
	}
}