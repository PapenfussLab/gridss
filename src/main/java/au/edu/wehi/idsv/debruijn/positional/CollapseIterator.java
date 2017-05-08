package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.Set;
import java.util.TreeSet;

import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Sets;

import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNodeUtil;
import it.unimi.dsi.fastutil.ints.IntRBTreeSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;

/**
 * Base class for path collapse logic
 * 
 * Records are processed within the processing window
 *    + post margin                                               + pre margin
 *     |----maxPathCollapseLength----|----maxPathCollapseLength----|
 *    |-------------------------------------------------------------|
 *    ^                                                             ^
 *   emit                                                          input
 * position                                                       position
 * 
 * 
 * @author Daniel Cameron
 *
 */
public abstract class CollapseIterator implements PeekingIterator<KmerPathNode> {
	private final PeekingIterator<KmerPathNode> underlying;
	protected final int k;
	private final int maxCollapseLength;
	protected final int maxBasesMismatch;
	private final NavigableSet<KmerPathNode> processed = new TreeSet<KmerPathNode>(KmerNodeUtil.ByFirstStartEndKmerReference);
	private final NavigableSet<KmerPathNode> unprocessed = new TreeSet<KmerPathNode>(KmerNodeUtil.ByLastEndStartKmerReference);
	private final int processOffset;
	private int lastEmitPosition = Integer.MIN_VALUE / 2; // moved away from MIN_VALUE to prevent underflow the calculating collapse window size
	protected int inputPosition = Integer.MIN_VALUE / 2;
	protected long nodesTraversed = 0;
	protected long leavesCollapsed = 0;
	protected long branchesCollapsed = 0;
	protected long consumed =  0;
	protected int maxNodeWidth = 0;
	protected int maxNodeLength = 0;
	protected int postCollapseBufferSize;
	private int emitOffset() {
		// records ending before this position cannot be changed by subsequent operations
		int unchangedOffset = processOffset + maxNodeLength + maxNodeWidth + maxCollapseLength + 1
				+ postCollapseBufferSize; // add extra buffer to allow for chained collapsing
		// records are ready to emit when their last kmer end position position is before unchangedOffset
		// we need to resort so they are output in the correct start kmer order
		return unchangedOffset + maxNodeLength + maxNodeWidth + 1;
	}
	protected int currentProcessPosition() {
		return inputPosition - processOffset;
	}
	protected int currentEmitPosition() {
		return inputPosition - emitOffset();
	}
	/**
	 * Indicate whether merged nodes should be further considered for collapse
	 * 
	 * Note: a leaf with greater weight that could be merged into a merged
	 * path will still not be reconsidered since the first leaf node is
	 * not actually on path.
	 * @return
	 */
	protected abstract boolean reprocessMergedNodes();
	public CollapseIterator(
			Iterator<KmerPathNode> it,
			int k,
			int maxPathCollapseLength,
			int maxBasesMismatch,
			int preCollapseBufferSize,
			int postCollapseBufferSize) {
		this.underlying = Iterators.peekingIterator(it);
		this.k = k;
		this.maxBasesMismatch = maxBasesMismatch;
		this.maxCollapseLength = maxPathCollapseLength;
		this.processOffset = maxPathCollapseLength + 1 + preCollapseBufferSize;
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
		lastEmitPosition = node.firstStart();
		return node;
	}
	@Override
	public KmerPathNode peek() {
		ensureBuffer();
		KmerPathNode node = processed.first();
		lastEmitPosition = node.firstStart();
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
			while (collapse()) {
				// collapse as much as we can
			} 
		}
	}
	/**
	 * Adds all nodes with first kmer starting at or before the current inputPosition
	 * to processing queue
	 */
	private void loadGraphNodes() {
		while (underlying.hasNext() && underlying.peek().firstStart() <= inputPosition) {
			KmerPathNode node = underlying.next();
			consumed++;
			maxNodeWidth = Math.max(maxNodeWidth, node.width());
			maxNodeLength = Math.max(maxNodeLength, node.length());
			unprocessed.add(node);
		}
	}
	private boolean collapse() {
		int collapseCount = 0;
		while (!unprocessed.isEmpty() && unprocessed.first().lastEnd() < currentProcessPosition()) {
			// calculate how much we can collapse
			int emitCollapseMargin = unprocessed.first().firstStart() - lastEmitPosition - 1;
			int inputCollapseMargin = inputPosition - unprocessed.first().lastEnd() - 1;
			int collapseWidth = Math.min(Math.min(emitCollapseMargin, inputCollapseMargin), maxCollapseLength);
			if (collapseNext(collapseWidth)) {
				collapseCount++;
			}
		}
		return collapseCount > 0;
	}
	private boolean collapseNext(int maxCollapseLength) {
		KmerPathNode node = unprocessed.pollFirst();
		processed.add(node);
		return collapse(node, maxCollapseLength);
	}
	/**
	 * Collapses paths involving the given node
	 * @param node
	 * @return
	 */
	protected abstract boolean collapse(KmerPathNode node, int maxCollapseLength);
	/**
	 * Merges the given source path into the target path 
	 * @param sourcePath path to merge
	 * @param targetPath merge destination
	 * @param sourceSkipKmers number of starting kmers to ignore in source
	 * @param targetSkipKmers number of starting kmers to ignore in target
	 */
	protected void merge(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath, int sourceSkipKmers, int targetSkipKmers) {
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
	/**
	 * Trims common path prefix and suffix nodes from both paths 
	 * @param path1
	 * @param path2
	 */
	private void trimCommon(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath) {
		while (!sourcePath.isEmpty() && !targetPath.isEmpty() && sourcePath.get(0).equals(targetPath.get(0))) {
			sourcePath.remove(0);
			targetPath.remove(0);
		}
		while (!sourcePath.isEmpty() && !targetPath.isEmpty() && sourcePath.get(sourcePath.size() - 1).equals(targetPath.get(targetPath.size() - 1))) {
			sourcePath.remove(sourcePath.size() - 1);
			targetPath.remove(targetPath.size() - 1);
		}
	}
	private void merge(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath) {
		assert(sourcePath.get(0).width() == targetPath.get(0).width());
		assert(sourcePath.get(0).firstStart() == targetPath.get(0).firstStart());
		assert(DeBruijnSequenceGraphNodeUtil.basesDifferent(k, sourcePath, targetPath) <= maxBasesMismatch);
		trimCommon(sourcePath, targetPath);
		assertKmerPathNodesUnique(sourcePath, targetPath);
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
		assert(source.size() <= target.size());
		// merge the common nodes
		for (int i = 0; i < source.size(); i++) {
			KmerPathNode toMerge = source.get(i);
			KmerPathNode into = target.get(i);
			merge(toMerge, into);
		}
		if (reprocessMergedNodes()) {
			for (int i = 0; i < source.size(); i++) {
				processed.remove(target.get(i));
				unprocessed.add(target.get(i));
			}
		}
	}
	private void assertKmerPathNodesUnique(List<KmerPathSubnode> sourcePath, List<KmerPathSubnode> targetPath) {
		Set<KmerPathNode> unique = Sets.newIdentityHashSet();
		for (KmerPathSubnode sn : Iterables.concat(sourcePath, targetPath)) {
			unique.add(sn.node());
		}
		assert(unique.size() == sourcePath.size() + targetPath.size());
	}
	/**
	 * Ensures that the 
	 * @param splitAfter splits after the given number of bases
	 * @param path path to split
	 */
	private void lengthSplit(int splitAfter, List<KmerPathSubnode> path) {
		assert(splitAfter >= 0);
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
				// add the split immediately before our node
				path.add(index, new KmerPathSubnode(split, n.firstStart(), n.firstEnd()));
				// and update the bounds on our node to their new starting position since the start position has shifted
				path.set(index + 1, new KmerPathSubnode(n.node(), n.firstStart() + splitLength, n.firstEnd() + splitLength));
				return;
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
		if (toMerge == into) return; // nothing to do
		assert(toMerge.lastStart() == into.lastStart());
		assert(toMerge.lastEnd() == into.lastEnd());
		assert(toMerge.length() == into.length());
		processed.remove(toMerge);
		unprocessed.remove(toMerge);
		// merge nodes
		into.merge(toMerge);
	}
	public int tracking_processedSize() {
		return processed.size();
	}
	public int tracking_inputPosition() {
		return inputPosition;
	}
	public long tracking_underlyingConsumed() {
		return consumed;
	}
	public int tracking_unprocessedSize() {
		return unprocessed.size();
	}
	public long tracking_traversalCount() {
		return nodesTraversed;
	}
	public long tracking_branchCollapseCount() {
		return branchesCollapsed;
	}
	public long tracking_leafCollapseCount() {
		return leavesCollapsed;
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}