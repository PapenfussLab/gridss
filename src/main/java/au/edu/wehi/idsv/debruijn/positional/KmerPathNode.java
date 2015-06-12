package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class KmerPathNode implements KmerNode, DeBruijnSequenceGraphNode {
	private static final LongArrayList EMPTY_KMER_LIST = new LongArrayList();
	private static final List<KmerPathNode> EMPTY_EDGE_LIST = ImmutableList.of();
	private LongArrayList kmers; // FIXME: replace with 2-bit encoding of kmer sequence
	private LongArrayList additionalKmers = null;
	private IntArrayList weight;
	private int totalWeight;
	private int start;
	private int end;
	private boolean reference;
	private int versionId = 0;
	/**
	 * Edges
	 */
	private ArrayList<KmerPathNode> nextList = null;
	private ArrayList<KmerPathNode> prevList = null;
	private boolean edgesSorted = true;
	/**
	 * Final kmer in path graph
	 */
	public long kmer() { return kmer(length() - 1); }
	/**
	 * First possible position of final kmer
	 */
	public int startPosition() { return startPosition(length() - 1); }
	/**
	 * Last possible position of final kmer
	 */
	public int endPosition() { return endPosition(length() - 1); }
	public long kmer(int offset) { return kmers.get(offset); }
	public int startPosition(int offset) { return start + offset; }
	public int endPosition(int offset) { return end + offset; }
	public int weight() { return totalWeight; }
	@Override
	public int weight(int offset) {
		return weight.getInt(offset);
	}
	public boolean isReference() { return reference; }
	public int length() { return kmers.size(); }
	/**
	 * Structural version identifier. The version number is changed
	 * whenever a structural modification is made to the path node.
	 * 
	 * Non-structural modifications such as weight change do not change the version.
	 * 
	 * Changes to edges are not considered structural changes to the node itself. 
	 * 
	 * @return incrementing version number
	 */
	public int version() { return versionId; }
	/**
	 * List of kmers that have been collapsed into this path
	 * 
	 * Note: this list is unordered
	 * @return collapsed kmers. Callers must not modify the returned list
	 */
	public LongArrayList collapsedKmers()
	{
		return additionalKmers != null ? additionalKmers : EMPTY_KMER_LIST;
	}
	public KmerPathNode(long kmer, int weight, int start, int end, boolean reference) {
		this.kmers = new LongArrayList(1);
		this.kmers.add(kmer);
		this.weight = new IntArrayList(1);
		this.weight.add(weight);
		this.totalWeight = weight;
		this.start = start;
		this.end = end;
		this.reference = reference;
	}
	private KmerPathNode(LongArrayList kmer, IntArrayList weight, int totalWeight, int start, int end, boolean reference) {
		this.kmers = kmer.clone();
		this.weight = weight.clone();
		this.totalWeight = totalWeight;
		this.start = start;
		this.end = end;
		this.reference = reference;
	}
	private KmerPathNode(LongArrayList kmer, IntArrayList weight, int start, int end, boolean reference) {
		this(kmer, weight, sumWeights(weight), start, end, reference);
	}
	public KmerPathNode(KmerNode node) {
		this(node.kmer(), node.weight(), node.startPosition(), node.endPosition(), node.isReference());
	}
	private static int sumWeights(IntArrayList weight) {
		int sum = 0;
		for (int i = 0; i < weight.size(); i++) {
			sum += weight.getInt(i);
		}
		return sum;
	}
	public void append(KmerAggregateNode node) {
		assert(node.startPosition() == startPosition(length() - 1) + 1);
		assert(node.endPosition() == endPosition(length() - 1) + 1);
		assert(node.isReference() == isReference());
		kmers.add(node.kmer());
		weight.add(node.weight());
		totalWeight += node.weight();
		reference |= node.isReference();
		versionId++;
	}
	/**
	 * Adds the given node to the front of this one
	 * @param node
	 */
	public void prepend(KmerPathNode node) {
		assert(startPosition() == node.startPosition(length() - 1) + 1);
		assert(endPosition() == node.endPosition(length() - 1) + 1);
		assert(isReference() == node.isReference());
		assert(node.next().size() == 1);
		assert(node.next().get(0) == this);
		assert(prev().size() == 1);
		assert(prev().get(0) == node);
		node.kmers.addAll(kmers);
		kmers = node.kmers;
		node.weight.addAll(weight);
		weight = node.weight;
		totalWeight += node.totalWeight;
		reference |= node.reference;
		if (node.additionalKmers != null) {
			if (additionalKmers == null) {
				additionalKmers = new LongArrayList(node.additionalKmers.size());
			}
			additionalKmers.addAll(node.additionalKmers);
		}
		nextList = node.nextList;
		versionId++;
		node.invalidate();
	}
	/**
	 * Merges the given node into this one
	 * @param toMerge
	 */
	public void merge(KmerPathNode toMerge) {
		assert(toMerge.startPosition() == startPosition());
		assert(toMerge.endPosition() == endPosition());
		assert(toMerge.length() == length());
		reference |= toMerge.reference;
		if (additionalKmers == null) {
			additionalKmers = new LongArrayList(toMerge.kmers.size() + toMerge.collapsedKmers().size());
		}
		additionalKmers.addAll(toMerge.kmers);
		additionalKmers.addAll(toMerge.collapsedKmers());
		totalWeight += toMerge.totalWeight;
		for (int i = 0; i < weight.size(); i++) {
			weight.set(i, weight.get(i) + toMerge.weight.get(i));
		}
		versionId++;
		replaceEdges(toMerge, this);
		toMerge.invalidate();
	}
	/**
	 * Indicates that this node has been transformed or deleted
	 *  
	 *  Edges 
	 */
	public void invalidate() {
		kmers = null;
		weight = null;
		nextList = null;
		prevList = null;
		additionalKmers = null;
		versionId++;
	}
	public boolean isValid() {
		return kmers != null;
	}
	/**
	 * Successor nodes, ordered by adjacency position  
	 * 
	 */
	public List<KmerPathNode> next() {
		if (nextList == null) return EMPTY_EDGE_LIST;
		ensureEdgesSorted();
		return nextList;
	}
	/**
	 * Precedessor nodes, ordered by adjacency position  
	 * 
	 */
	public List<KmerPathNode> prev() {
		if (prevList == null) return EMPTY_EDGE_LIST;
		ensureEdgesSorted();
		return prevList;
	}
	public static void addEdge(KmerPathNode from, KmerPathNode to) {
		assert(from.isValid());
		assert(to.isValid());
		assert(!from.next().contains(to));
		assert(!to.prev().contains(from));
		from.nextList.add(to);
		to.prevList.add(from);
		if (from.nextList.size() > 1 && KmerPathNode.ByFirstKmerStartPosition.compare(from.nextList.get(from.nextList.size() - 2), from.nextList.get(from.nextList.size() - 1)) > 0) {
			from.edgesSorted = false;
		}
		if (to.prevList.size() > 1 && KmerPathNode.ByFirstKmerStartPosition.compare(to.prevList.get(to.prevList.size() - 2), to.prevList.get(to.prevList.size() - 1)) > 0) {
			to.edgesSorted = false;
		}
	}
	/**
	 * Replaces all edges coming to/from the source node, with edges to/from the target node
	 * @param toRemoveFrom source node to remove all edges from
	 * @param toAddTo target node to replace with
	 */
	private static void replaceEdges(KmerPathNode toRemoveFrom, KmerPathNode toAddTo) {
		for (KmerPathNode n : toRemoveFrom.next()) {
			if (!toAddTo.next().contains(n)) {
				addEdge(toAddTo, n);
			}
		}
		for (KmerPathNode n : toRemoveFrom.prev()) {
			if (!toAddTo.prev().contains(n)) {
				addEdge(n, toAddTo);
			}
		}
	}
	private void ensureEdgesSorted() {
		if (!edgesSorted) {
			if (nextList != null) {
				Collections.sort(nextList, KmerPathNode.ByFirstKmerStartPosition);
			}
			if (prevList != null) {
				Collections.sort(prevList, KmerPathNode.ByStartPosition);
			}
			edgesSorted = true;
		}
	}
	/**
	 * Splits out a the node by creating a new node containing the first given number of bases
	 * 
	 * Note: PositionalDeBruijnGraphPathNodeGraphSimplifier relies on the fact that the node
	 * is truncated at the start, not the end.
	 * 
	 * @param length length of node to split out
	 * @return new predecessor node
	 */
	public KmerPathNode splitAtLength(int firstNodeLength) {
		if (firstNodeLength == 0 || firstNodeLength == length()) return this;
		assert(firstNodeLength > 0);
		assert(firstNodeLength < length());
		// copy our new kmers and weights
		LongArrayList kmerSecond = new LongArrayList(kmers.subList(firstNodeLength, length()));
		IntArrayList weightSecond = new IntArrayList(weight.subList(firstNodeLength, length()));
		// let split own our current arrays
		this.kmers.removeElements(firstNodeLength, this.kmers.size());
		this.weight.removeElements(firstNodeLength, this.weight.size());
		KmerPathNode split = new KmerPathNode(
				this.kmers,
				this.weight,
				start,
				end,
				reference);
		// Update incoming edges to split
		ArrayList<KmerPathNode> prevEmpty = split.prevList;
		split.prevList = this.prevList;
		this.prevList = prevEmpty;
		for (KmerPathNode n : split.prevList) {
			replaceFirst(n.nextList, this, split);
		}
		addEdge(split, this);
		// edge sorting invariant remains unchanged
				
		this.kmers = kmerSecond;
		this.weight = weightSecond;
		this.totalWeight -= split.weight();
		this.start += firstNodeLength;
		this.end += firstNodeLength;
		this.versionId++;
		return split;
	}
	private static void replaceFirst(List<KmerPathNode> list, KmerPathNode existing, KmerPathNode replaceWith) {
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i) == existing) { 
				list.set(i, replaceWith);
				return;
			}
		}
		throw new IllegalStateException("Could not replace non-existant element " + existing.toString());
	}
	/**
	 * Splits at the given internal start position
	 * 
	 * @param internalStart start position to split at
	 * @return new node ending immediately before the given position
	 */
	public KmerPathNode splitAtStartPosition(int newStartPosition) {
		assert(newStartPosition > start);
		assert(newStartPosition < end);
		KmerPathNode split = new KmerPathNode(kmers, weight, totalWeight, start, newStartPosition - 1, reference);
		this.start = newStartPosition;
		
		ArrayList<KmerPathNode> newNextThis = new ArrayList<KmerPathNode>(nextList.size());
		ArrayList<KmerPathNode> newNextSplit = new ArrayList<KmerPathNode>(nextList.size());
		for (KmerPathNode adj : nextList) {
			if (IntervalUtil.overlapsClosed(adj.startPosition(0) - 1, adj.endPosition(0) - 1, this.start, this.end)) {
				newNextThis.add(adj);
			}
			if (IntervalUtil.overlapsClosed(adj.startPosition(0) - 1, adj.endPosition(0) - 1, split.start, split.end)) {
				newNextSplit.add(adj);
			}
		}
		this.nextList = newNextThis;
		split.nextList = newNextSplit;
		
		ArrayList<KmerPathNode> newPrevThis = new ArrayList<KmerPathNode>(prevList.size());
		ArrayList<KmerPathNode> newPrevSplit = new ArrayList<KmerPathNode>(prevList.size());
		for (KmerPathNode adj : prevList) {
			if (IntervalUtil.overlapsClosed(adj.startPosition(0) + 1, adj.endPosition(0) + 1, this.start, this.end)) {
				newPrevThis.add(adj);
			}
			if (IntervalUtil.overlapsClosed(adj.startPosition(0) + 1, adj.endPosition(0) + 1, split.start, split.end)) {
				newPrevSplit.add(adj);
			}
		}
		this.prevList = newPrevThis;
		split.prevList = newPrevSplit;
		split.edgesSorted = this.edgesSorted; // edge list order retained
		this.versionId++;
		return split;
	}
	public boolean sanityCheck(int k, int maxSupportWidth, int maxPathLength) {
		assert(isValid());
		assert(start <= end);
		assert(totalWeight > 0);
		assert(length() <= maxPathLength);
		assert(end - start <= maxSupportWidth);
		assert(kmers.size() == length());
		assert(weight.size() == length());
		for (int i = 1; i < length(); i++) {
			assert(KmerEncodingHelper.isNext(k, kmers.getLong(i - 1), kmers.getLong(i)));
		}
		assert(sumWeights(weight) == totalWeight);
		for (KmerPathNode next : nextList) {
			assert(KmerEncodingHelper.isNext(k, kmer(), next.kmer(0)));
			assert(IntervalUtil.overlapsClosed(startPosition() + 1, endPosition() + 1, next.startPosition(0), next.endPosition(0)));
			assert(next.isValid());
		}
		for (KmerPathNode prev : prevList) {
			assert(KmerEncodingHelper.isNext(k, prev.kmer(), kmer(0)));
			assert(IntervalUtil.overlapsClosed(startPosition(0) - 1, endPosition(0) - 1, prev.startPosition(), prev.endPosition()));
			assert(prev.isValid());
		}
		if (edgesSorted) {
			for (int i = 1; i < nextList.size(); i++) {
				assert(nextList.get(i - 1).startPosition(0) < nextList.get(i).startPosition(0));
			}
			for (int i = 1; i < prevList.size(); i++) {
				assert(prevList.get(i - 1).startPosition() < prevList.get(i).startPosition());
			}
		}
		assert(EMPTY_KMER_LIST.size() == 0); // fastutil doesn't have ImmutableList wrappers
		return true;
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s%d, %d kmers", isReference() ? "R" : " ",  startPosition(), endPosition(), weight(), length());
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + kmers.get(0).hashCode();
		result = prime * result + start;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		KmerPathNode other = (KmerPathNode) obj;
		if (end != other.end)
			return false;
		if (kmers.get(0) != other.kmers.get(0))
			return false;
		if (start != other.start)
			return false;
		return true;
	}
	public static final Ordering<KmerPathNode> ByFirstKmerStartPosition = new Ordering<KmerPathNode>() {
		@Override
		public int compare(KmerPathNode left, KmerPathNode right) {
			return Ints.compare(left.startPosition(0), right.startPosition(0));
		}
	};
	public static final Ordering<KmerPathNode> ByFirstKmerEndPosition = new Ordering<KmerPathNode>() {
		@Override
		public int compare(KmerPathNode left, KmerPathNode right) {
			return Ints.compare(left.startPosition(0), right.startPosition(0));
		}
	};
}
