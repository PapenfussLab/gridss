package au.edu.wehi.idsv.debruijn.positional;

import htsjdk.samtools.util.RuntimeEOFException;
import it.unimi.dsi.fastutil.ints.IntRBTreeSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.Stack;
import java.util.TreeSet;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNodeUtil;
import au.edu.wehi.idsv.debruijn.HighestWeightSimilarPath;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
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
 * For simplicity, each node processed for all collapse types once:
 *                             |----maxSupportWidth+maxPathLength----|
 *    |----maxcollaspeLength---|                                     |----maxcollaspeLength---|
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
	private final AssemblyParameters ap;
	private final int maxSupportWidth;
	private final int processOffset;
	private final int emitOffset;
	private final int maxCollapseLength;
	private final NavigableSet<KmerPathNode> processed = new TreeSet<KmerPathNode>(KmerPathNode.ByFirstKmerStartPosition);
	private final NavigableSet<KmerPathNode> unprocessed = new TreeSet<KmerPathNode>(KmerPathNode.ByEndPositionKmer);
	private int inputPosition = Integer.MIN_VALUE;
	public PathCollapseIterator(
			Iterator<KmerPathNode> it,
			AssemblyParameters ap,
			int maxSupportWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.ap = ap;
		this.maxSupportWidth = maxSupportWidth;
		int maxPathLength = ap.positionalMaxPathLengthInBases(maxSupportWidth);
		this.maxCollapseLength = ap.positionalMaxPathCollapseLengthInBases(maxSupportWidth);
		this.processOffset = maxCollapseLength + 1;
		// records before this position cannot be changed by subsequent operations
		int unchangedOffset = processOffset + maxPathLength + maxSupportWidth + maxCollapseLength + 1;
		// records are added to the emit queue when their last kmer end position position is before unchangedOffset
		// since we need to resort so they are output in the correct order
		this.emitOffset = unchangedOffset - maxSupportWidth; 
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !processed.isEmpty();
	}
	@Override
	public KmerPathNode next() {
		ensureBuffer();
		KmerPathNode node = processed.pollFirst();
		return node;
	}
	private void ensureBuffer() {
		while (inputPosition < Integer.MAX_VALUE && (processed.isEmpty() || processed.first().startPosition(0) < inputPosition - emitOffset)) {
			// advance graph position
			if (underlying.hasNext()) {
				inputPosition = underlying.peek().startPosition(0);
				loadGraphNodes();
			} else {
				inputPosition = Integer.MAX_VALUE;
			}
			while (collapse() > 0) { } // collapse as much as we can
		}
	}
	/**
	 * Adds all nodes with first kmer starting at or before the current inputPosition
	 * to processing queue
	 */
	private void loadGraphNodes() {
		while (underlying.hasNext() && underlying.peek().startPosition(0) <= inputPosition) {
			unprocessed.add(underlying.next());
		}
	}
	private int collapse() {
		int collapseCount = 0;
		while (!unprocessed.isEmpty() && unprocessed.first().endPosition() < inputPosition - processOffset) {
			collapseCount += collapseNext();
		}
		return collapseCount;
	}
	private int collapseNext() {
		KmerPathNode node = unprocessed.pollFirst();
		processed.add(node);
		return collapse(node);
	}
	private int collapse(KmerPathNode node) {
		KmerPathSubnode root = new KmerPathSubnode(node); 
		// TODO: traverse backward looking for:
		// leaves to collapse
		// paths to collapse
		List<KmerPathSubnode> nextNodes = root.next();
		for (int i = 0; i < nextNodes.size(); i++) {
			for (int j = 1; j < nextNodes.size(); j++) {
				// Alternate algorithm:
			}
		}
		// traverse forward looking for:
		// leaves to collapse
		// paths already covered by backward traverse from the other end
		throw new RuntimeException("NYI");
	}
	private static List<List<KmerPathSubnode>> findSimilarPaths(KmerPathNode root, boolean findLeaf, boolean findCommonChild, int maxLength, boolean traverseForward) {
		// node is on leaf path if <= 1 prev and next edge (and previous node was on leaf path) 
		
		// TODO: dynamic programming algorithm:
		// Common child:
		// find all reachable nodes in both graphs
		// find overlaps
		
		// if too different return
		// if terminal leaf return
		// if have matching end node, return
		
		// while path queues not empty
			// advance the shorter of the two paths
		throw new RuntimeException("NYI");
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
			KmerPathNode split = lengthSplit(current, startPosition - current.startPosition(0));
			result.add(split);
			assert(current.startPosition(0) == startPosition);
		}
		result.add(current);
		while (it.hasNext()) {
			result.add(it.next());
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
		int preLength = n.firstKmerStartPosition() - pn.startPosition(0);
		if (preLength != 0) {
			NavigableSet<KmerPathNode> queue = processed.contains(pn) ? processed : unprocessed;
			queue.remove(pn);
			KmerPathNode preNode = pn.splitAtLength(preLength);
			queue.add(pn);
			queue.add(preNode);
		}
		int postLength = pn.endPosition(0) - n.firstKmerEndPosition();
		if (postLength != 0) {
			NavigableSet<KmerPathNode> queue = processed.contains(pn) ? processed : unprocessed;
			queue.remove(pn);
			KmerPathNode midNode = pn.splitAtLength(preLength);
			queue.add(pn);
			queue.add(midNode);
			pn = midNode;
		}
		return pn;
	}
	private void merge(KmerPathNode toMerge, KmerPathNode into) {
		assert(toMerge.startPosition() == into.startPosition());
		assert(toMerge.endPosition() == into.endPosition());
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