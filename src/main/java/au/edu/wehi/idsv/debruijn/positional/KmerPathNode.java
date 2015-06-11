package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class KmerPathNode extends KmerNode implements DeBruijnSequenceGraphNode {
	private LongArrayList kmers; // FIXME: replace with 2-bit encoding of kmer sequence
	private IntArrayList weight;
	private int totalWeight;
	private int start;
	private int end;
	private boolean reference;
	private int versionId = 0;
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
	 * Indicates that this node has been transformed or deleted
	 *  
	 *  Edges 
	 */
	public void invalidate() {
		kmers = null;
		weight = null;
		versionId++;
	}
	public boolean isValid() {
		return kmers != null;
	}
	/**
	 * Edges
	 */
	private ArrayList<KmerPathNode> nextList = new ArrayList<KmerPathNode>(3);
	private ArrayList<KmerPathNode> prevList = new ArrayList<KmerPathNode>(3);
	/**
	 * Successor nodes, ordered by adjacency position  
	 * 
	 */
	public List<KmerPathNode> next() {
		ensureEdgesSorted();
		return nextList;
	}
	/**
	 * Precedessor nodes, ordered by adjacency position  
	 * 
	 */
	public List<KmerPathNode> prev() {
		ensureEdgesSorted();
		return prevList;
	}
	private boolean edgesSorted = false;
	public static void addEdge(KmerPathNode from, KmerPathNode to) {
		assert(from.isValid());
		assert(to.isValid());
		from.nextList.add(to);
		to.prevList.add(from);
		from.edgesSorted = false;
		to.edgesSorted = false;
	}
	private void ensureEdgesSorted() {
		if (!edgesSorted) {
			Collections.sort(nextList, KmerPathNode.ByFirstKmerStartPosition);
			Collections.sort(prevList, KmerPathNode.ByStartPosition);
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
