package au.edu.wehi.idsv.debruijn.positional;

import java.util.HashMap;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;

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
 * width of the merge is bounded by maxSupportWidth + maxPathLength = maxNodeSize in either direction
 * 
 * @author cameron.d
 *
 */
public class PathSimplificationIterator implements Iterator<KmerPathNode> {
	private final HashMap<KmerNode.KmerStartPositionEqualityWrapper, KmerPathNode> startLookup = new HashMap<KmerNode.KmerStartPositionEqualityWrapper, KmerPathNode>();
	private final HashMap<KmerNode.KmerEndPositionEqualityWrapper, KmerPathNode> endLookup = new HashMap<KmerNode.KmerEndPositionEqualityWrapper, KmerPathNode>();
	/**
	 *  |---- maxNodeLength+maxNodeWidth----|                                       
	 *                                      |---- maxNodeLength+maxNodeWidth----|
	 *  ^                                   ^                                   ^           
	 *  |                                  process                         input position
	 * unchanged                          position
	 *  offset
	 *  
	 *  Worst case scenario results in merge of node of
	 *  kmer length of 1, support width of 1 until 
	 *  node is maxNodeLength + maxNodeWidth in size 
	 * 
	 */
	private final NavigableSet<KmerPathNode> unprocessed = new TreeSet<KmerPathNode>(KmerPathNode.ByFirstKmerStartPosition);
	private final NavigableSet<KmerPathNode> processed = new TreeSet<KmerPathNode>(KmerPathNode.ByFirstKmerStartPosition);
	private final int maxNodeLength;
	private final int maxNodeWidth;
	private final int processOffset;
	private final int emitOffset;
	private final PeekingIterator<KmerPathNode> underlying;
	private int inputPosition;
	public PathSimplificationIterator(
			Iterator<KmerPathNode> it,
			int maxNodeLength,
			int maxNodeWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxNodeLength = maxNodeLength;
		this.maxNodeWidth = maxNodeWidth;
		this.processOffset = maxNodeLength + maxNodeWidth + 1;
		this.emitOffset = processOffset + maxNodeLength + maxNodeWidth + 1;
	}
	private void enqueueForProcessing(KmerPathNode node) {
		assert(node.isValid());
		startLookup.put(new KmerNode.KmerStartPositionEqualityWrapper(node), node);
		endLookup.put(new KmerNode.KmerEndPositionEqualityWrapper(node), node);
		unprocessed.add(node);
	}
	private void dequeue(KmerPathNode node) {
		startLookup.remove(new KmerNode.KmerStartPositionEqualityWrapper(node));
		endLookup.remove(new KmerNode.KmerEndPositionEqualityWrapper(node));
		unprocessed.remove(node);
		processed.remove(node);
	}
	private boolean compress(KmerPathNode node) {
		// merge with successor
		if (mergeWithNextKmers(node)) return true;
		if (node.prev().size() == 1) {
			// merge with predecessor
			if (mergeWithNextKmers(node.prev().get(0))) return true;
		}
		if (mergeAdjacent(endLookup.get(new KmerNode.KmerEndPositionEqualityLookup(node.kmer(), node.startPosition() - 1)), node)) return true;
		if (mergeAdjacent(node, startLookup.get(new KmerNode.KmerStartPositionEqualityLookup(node.kmer(), node.endPosition() + 1)))) return true;
		return false;
	}
	private boolean mergeAdjacent(KmerPathNode first, KmerPathNode second) {
		if (first == null || second == null) return false;
		assert(first.startPosition() <= second.startPosition());
		if (first.canCoalese(second)
				&& first.width() + second.width() <= maxNodeWidth) {
			dequeue(first);
			dequeue(second);
			second.coaleseAdjacent(first);
			enqueueForProcessing(second);
			return true;
		}
		return false;
	}
	private boolean mergeWithNextKmers(KmerPathNode node) {
		if (node.next().size() == 1) {
			KmerPathNode candidate = node.next().get(0);
			if (candidate.prev().size() == 1
					&& candidate.startPosition(0) == node.endPosition() + 1
					&& candidate.endPosition(0) == node.endPosition() + 1
					&& node.length() + candidate.length() <= maxNodeLength) {
				// we can concatenate these nodes together
				dequeue(node);
				dequeue(candidate);
				candidate.prepend(node);
				enqueueForProcessing(candidate);
				return true;
			}
		}
		return false;
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
		assert(node.startPosition(0) < inputPosition - emitOffset);
		dequeue(node);
		return node;
	}
	private void ensureBuffer() {
		while (inputPosition < Integer.MAX_VALUE && (processed.isEmpty() || processed.first().startPosition(0) > inputPosition - emitOffset)) {
			// advance graph position
			if (underlying.hasNext()) {
				inputPosition = underlying.peek().startPosition(0);
			} else {
				inputPosition = Integer.MAX_VALUE;
			}
			loadGraphNodes();
			process();
		}
	}
	private void process() {
		while (!unprocessed.isEmpty() && unprocessed.first().startPosition(0) <= inputPosition - processOffset) {
			KmerPathNode node = unprocessed.first();
			if (!compress(node)) {
				unprocessed.remove(node);
				processed.add(node);
			}
		}
	}
	private void loadGraphNodes() {
		while (underlying.hasNext() && underlying.peek().startPosition(0) <= inputPosition) {
			enqueueForProcessing(underlying.next());
		}
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
