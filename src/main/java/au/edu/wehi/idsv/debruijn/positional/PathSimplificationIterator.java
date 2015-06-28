package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.objects.ObjectOpenCustomHashSet;

import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.TreeSet;

import au.edu.wehi.idsv.Defaults;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Simplifies the graph until the number of nodes used to represent the graph cannot be further reduced
 * 
 *  Without loss of generality, consider merging two nodes N_1, N_2 with N_1 \in next(M): A holds before merge
 *  - merging nodes N_1 and N_2 into N' will cause A to be violated iff:
 *   - adjacent node M can be collapsed into N', but could not be collapsed into N_1.
 * To merge nodes:
 *  start(N_1) = start(N_2)
 *  end(N_1) = end(N_2)
 *  kmers with hamming distance threshold
 *  |next(N_12)| >= |next(N_1)| degree does not decrease
 *  M degree reduced by 1 iff N_2 \in next(M)
 *  	N_12 can be merged with M 
 *   *  
 *  Consider the following graph:
 *  [11-15] 1 GTAC (node A)
 *  [16-20] 2 GTAC (node B)
 *  [11-15] 1 GGAC (node C)
 *  [10-19] 2 GGTA (node D)
 *  
 *  We merges C into A:
 *  [11-15] 2 GTAC (node A,C)
 *  [16-20] 2 GTAC (node B)
 *  [10-19] 2 GGTA (node D)
 *  
 *  Merging AC with B (since they differ only in validity intervals, and these intervals are adjacent) results in:  
 *  [11-20] 2 GTAC (node A,B,C)
 *  [10-19] 2 GGTA (node D)
 *  
 *  Merging ABC with D:
 *  [10-19] 4 GGTAC (nodes A,B,C,D)
 *  
 * A larger example would result in additional nodes being further merged.
 * The number of possible merges is unbounded, but since each merge results in a node containing N_1,
 * width and length of the merge is bounded by maxSupportWidth and maxPathLength
 * 
 * @author cameron.d
 *
 */
public class PathSimplificationIterator implements Iterator<KmerPathNode> {
	private final ObjectOpenCustomHashSet<KmerPathNode> endLookup = new ObjectOpenCustomHashSet<KmerPathNode>(new KmerNodeUtil.HashByLastEndKmer<KmerPathNode>());
	/**
	 * Always merge with earlier nodes. Ordering by end position ensures that
	 * both kmer and position precedessors have been processed before each node
	 * is processed.
	 * 
	 *  |---- maxNodeLength+maxNodeWidth----|                                       
	 *                                      |
	 *  ^                                   ^           
	 *  |                                  input
	 * unchanged                          position
	 * offset
	 *  
	 *  Worst case scenario results in merge of node of
	 *  kmer length of 1, support width of 1 until 
	 *  node is maxNodeLength + maxNodeWidth in size
	 * 
	 */
	private final PriorityQueue<KmerPathNode> unprocessed = new PriorityQueue<KmerPathNode>(16, KmerNodeUtil.ByLastEnd);
	/**
	 * Nodes that have been processed, but could be modified further
	 * 
	 */
	private final TreeSet<KmerPathNode> processed = new TreeSet<KmerPathNode>(KmerNodeUtil.ByFirstStartKmer);
	private final int maxLength;
	private final int maxWidth;
	private final PeekingIterator<KmerPathNode> underlying;
	private int inputPosition;
	public PathSimplificationIterator(
			Iterator<KmerPathNode> it,
			int maxLength,
			int maxWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxLength = maxLength;
		this.maxWidth = maxWidth;
	}
	/**
	 * Determines whether it is possible for the given node to be expanded
	 * to contain a kmer at the given position
	 * @param node node to consider
	 * @param position position to check
	 * @return true if a merge is possible, false otherwise
	 */
	private boolean couldMergeToIncludeKmerAt(KmerPathNode node, int position) {
		int latestFirstKmerEnd = node.lastStart() + maxWidth - 1;
		int latestLastKmerEnd = latestFirstKmerEnd + maxLength - 1;
		return latestLastKmerEnd >= position;
	}
	/**
	 * Merges the given node with the a matching previous node.
	 * 
	 * @param node node to attempt merge
	 * @return true if the node was merged, false otherwise.
	 */
	private boolean reduce(KmerPathNode node) {
		KmerPathNode prev = prevKmerToMergeWith(node);
		if (prev != null) {
			assert(processed.contains(prev));
			processed.remove(prev);
			endLookup.remove(prev);
			node.prepend(prev);
			return true;
		}
		KmerPathNode adj = adjacentBeforeKmerToMergeWith(node);
		if (adj != null) {
			assert(processed.contains(adj));
			processed.remove(adj);
			endLookup.remove(adj);
			node.coaleseBeforeAdjacent(adj);
			return true;
		}
		return false;
	}
	private KmerPathNode adjacentBeforeKmerToMergeWith(KmerPathNode node) {
		KmerPathNode adj = endLookup.get(new KmerPathNode(node.lastKmer(), 0, node.lastStart() - 1, false, 0));
		if (adj != null
				&& node.canCoaleseBeforeAdjacent(adj)
				&& adj.width() + node.width() <= maxWidth) {
			return adj;
		}
		return null;
	}
	private KmerPathNode prevKmerToMergeWith(KmerPathNode node) {
		if (node.prev().size() == 1) {
			KmerPathNode prev = node.prev().get(0);
			if (prev.lastEnd() + 1 == node.firstKmerEnd()
					&& prev.lastStart() + 1 == node.firstKmerStart()
					&& prev.isReference() == node.isReference()
					&& prev.length() + node.length() <= maxLength
					&& prev.next().size() == 1) {
				return prev;
			}
		}
		return null;
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
		endLookup.remove(node);
		assert(!couldMergeToIncludeKmerAt(node, inputPosition));
		return node;
	}
	private void ensureBuffer() {
		while (inputPosition < Integer.MAX_VALUE && (processed.isEmpty() || couldMergeToIncludeKmerAt(processed.first(), inputPosition))) {
			// advance graph position
			if (underlying.hasNext()) {
				inputPosition = underlying.peek().firstKmerStart();
			} else {
				inputPosition = Integer.MAX_VALUE;
			}
			advance();
			process();
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheck());
			if (inputPosition == Integer.MAX_VALUE) {
				assert(!underlying.hasNext());
				assert(unprocessed.isEmpty());
			}
		}
	}
	private void process() {
		while (!unprocessed.isEmpty() && unprocessed.peek().lastEnd() <= inputPosition) {
			KmerPathNode node = unprocessed.poll();
			process(node);
		}
	}
	private void process(KmerPathNode node) {
		int beforeTotalWeight = 0;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			beforeTotalWeight = node.width() * node.weight() + processed.stream().mapToInt(n -> n.width() * n.weight()).sum();
			assert(!processed.contains(node));
		}
		while (reduce(node)) {
			// loop until we can't merge any more nodes into this one
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			int afterTotalWeight = node.width() * node.weight() + processed.stream().mapToInt(n -> n.width() * n.weight()).sum();
			assert(beforeTotalWeight == afterTotalWeight);
			assert(!processed.contains(node));
		}
		processed.add(node);
		endLookup.add(node);
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheck());
		}
	}
	/**
	 * Loads records from the underlying stream up to and including the current inputPosition.
	 */
	private void advance() {
		while (underlying.hasNext() && underlying.peek().firstKmerStart() <= inputPosition) {
			KmerPathNode nextRecord = underlying.next();
			assert(nextRecord.width() <= maxWidth);
			assert(nextRecord.length() <= maxLength);
			unprocessed.add(nextRecord);
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheck());
		}
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	public boolean sanityCheck() {
		assert(endLookup.size() == processed.size());
		for (KmerPathNode pn : processed) {
			// should not be able to reduce processed nodes any further
			assert(prevKmerToMergeWith(pn) == null);
			assert(adjacentBeforeKmerToMergeWith(pn) == null);
			assert(unprocessed.stream().filter(n -> pn.equals(n)).count() == 0);
		}
		for (KmerPathNode pn : unprocessed) {
			assert(processed.stream().filter(n -> pn.equals(n)).count() == 0);
		}
		return true;
	}
}
