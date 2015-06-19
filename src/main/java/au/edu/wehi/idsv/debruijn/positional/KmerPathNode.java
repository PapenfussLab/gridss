package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class KmerPathNode implements KmerNode, DeBruijnSequenceGraphNode {
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
	public static final Ordering<KmerPathNode> ByFirstKmerStartPositionKmer = new Ordering<KmerPathNode>() {
		@Override
		public int compare(KmerPathNode left, KmerPathNode right) {
			return ComparisonChain.start()
					.compare(left.startPosition(0), right.startPosition(0))
					.compare(left.kmer(), right.kmer())
					.result();
		}
	};	
	private static final LongArrayList EMPTY_KMER_LIST = new LongArrayList();
	private static final List<KmerPathNode> EMPTY_EDGE_LIST = ImmutableList.of();
	private static final Ordering<KmerPathNode> NEXT_SORT_ORDER = KmerPathNode.ByFirstKmerStartPosition;
	private static final Ordering<KmerNode> PREV_SORT_ORDER = KmerPathNode.ByStartPosition;
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
	public int width() { return end - start + 1; }
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
	public KmerPathNode(long kmer, int start, int end, boolean reference, int weight) {
		this.kmers = new LongArrayList(1);
		this.kmers.add(kmer);
		this.weight = new IntArrayList(1);
		this.weight.add(weight);
		this.totalWeight = weight;
		this.start = start;
		this.end = end;
		this.reference = reference;
	}
	private KmerPathNode(LongArrayList kmer, int start, int end, boolean reference, int totalWeight, IntArrayList weight) {
		this.kmers = kmer.clone();
		this.weight = weight.clone();
		this.totalWeight = totalWeight;
		this.start = start;
		this.end = end;
		this.reference = reference;
	}
	private KmerPathNode(LongArrayList kmer, int start, int end, boolean reference, IntArrayList weight) {
		this(kmer, start, end, reference, sumWeights(weight), weight);
	}
	public KmerPathNode(KmerNode node) {
		this(node.kmer(), node.startPosition(), node.endPosition(), node.isReference(), node.weight());
	}
	private static int sumWeights(IntArrayList weight) {
		int sum = 0;
		for (int i = 0; i < weight.size(); i++) {
			sum += weight.getInt(i);
		}
		return sum;
	}
	public void append(KmerNode node) {
		assert(!(node instanceof KmerPathNode)); // should be using prepend
		assert(node.startPosition() == startPosition() + 1);
		assert(node.endPosition() == endPosition() + 1);
		assert(node.isReference() == isReference());
		assert(nextList == null || nextList.size() == 0);
		//assert(KmerEncodingHelper.isNext(k, kmer(length() - 1), node.kmer())); 
		kmers.add(node.kmer());
		weight.add(node.weight());
		totalWeight += node.weight();
		reference |= node.isReference();
		versionId++;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			sanityCheck();
		}
	}
	/**
	 * Adds the given node to the front of this one
	 * @param node
	 */
	public void prepend(KmerPathNode node) {
		assert(startPosition(0) == node.startPosition() + 1);
		assert(endPosition(0) == node.endPosition() + 1);
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
		prevList = node.prevList;
		edgesSorted &= node.edgesSorted;
		if (prevList != null) {
			for (KmerPathNode n : prevList) {
				replaceFirst(n.nextList, node, this);
				n.edgesSorted = false;
			}
		}
		start = node.start;
		end = node.end;
		versionId++;
		node.prevList = null;
		node.nextList = null;
		node.invalidate();
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			sanityCheck();
		}
	}
	public boolean canCoaleseBeforeAdjacent(KmerPathNode node) {
		return start == node.end + 1
				&& length() == node.length()
				&& reference == node.reference
				&& totalWeight == node.totalWeight 
				&& kmers.equals(node.kmers)
				&& weight.equals(node.weight);
	}
	/**
	 * Merges the node covering the adjacent interval with matching kmers
	 * immediately preceeding this node in genomic 
	 * @param node 
	 */
	public void coaleseBeforeAdjacent(KmerPathNode node) {
		assert(canCoaleseBeforeAdjacent(node));
		replaceEdges(node, this);
		start = node.start;
		versionId++;
		node.invalidate();
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			sanityCheck();
		}
	}
	/**
	 * Merges the given alternate kmer path into this one
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
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			sanityCheck();
		}
	}
	/**
	 * Indicates that this node has been transformed or deleted
	 *  
	 *  Edges 
	 */
	public void invalidate() {
		assert(next().size() == 0);
		assert(prev().size() == 0);
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
	 * Divides this PathNode into Subnodes that each share the same PathNode sucessors
	 * @return KmerPathSubnodes fully covering this KmerPathNode
	 */
	public List<KmerPathSubnode> asSubnodesByNext() {
		List<KmerPathSubnode> subnodes = new ArrayList<KmerPathSubnode>(next().size() + 1);
		PriorityQueue<KmerPathNode> active = new PriorityQueue<KmerPathNode>(4, KmerPathNode.ByFirstKmerEndPosition);
		ensureEdgesSorted();
		if (nextList == null) {
			subnodes.add(new KmerPathSubnode(this));
		} else {
			int activeStart = start;
			int offset = 0;
			while (activeStart <= end) {
				while (offset < nextList.size() && nextList.get(offset).startPosition(0) - length() <= activeStart) {
					// add successors that should now be in scope
					active.add(nextList.get(offset++));
				}
				while (!active.isEmpty() && active.peek().endPosition(0) - length() < activeStart) {
					// remove out of scope successors
					active.poll();
				}
				int endPos = end;
				if (offset < nextList.size()) {
					int nextStart = nextList.get(offset).startPosition(0) - length();
					endPos = Math.min(endPos, nextStart - 1);
				}
				if (!active.isEmpty()) {
					int nextEnd = active.peek().endPosition(0) - length();
					endPos = Math.min(endPos, nextEnd);
				}
				subnodes.add(new KmerPathSubnode(this, activeStart, endPos));
				activeStart = endPos + 1;
			}
		}
		return subnodes;
	}
	/**
	 * Divides this PathNode into Subnodes that each share the same PathNode sucessors
	 * @return KmerPathSubnodes fully covering this KmerPathNode
	 */
	public List<KmerPathSubnode> asSubnodesByPrev() {
		List<KmerPathSubnode> subnodes = new ArrayList<KmerPathSubnode>(prev().size() + 1);
		PriorityQueue<KmerPathNode> active = new PriorityQueue<KmerPathNode>(4, KmerPathNode.ByFirstKmerEndPosition);
		ensureEdgesSorted();
		if (prevList == null) {
			subnodes.add(new KmerPathSubnode(this));
		} else {
			int activeStart = start;
			int offset = 0;
			while (activeStart <= end) {
				while (offset < prevList.size() && prevList.get(offset).start <= activeStart) {
					// add successors that should now be in scope
					active.add(prevList.get(offset++));
				}
				while (!active.isEmpty() && active.peek().endPosition() + 1 < activeStart) {
					// remove out of scope successors
					active.poll();
				}
				int endPos = end;
				if (offset < prevList.size()) {
					int nextStart = prevList.get(offset).startPosition() + 1;
					endPos = Math.min(endPos, nextStart - 1);
				}
				if (!active.isEmpty()) {
					int nextEnd = active.peek().endPosition() + 1;
					endPos = Math.min(endPos, nextEnd);
				}
				subnodes.add(new KmerPathSubnode(this, activeStart, endPos));
				activeStart = endPos + 1;
			}
		}
		return subnodes;
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
		if (from.nextList == null) {
			from.nextList = new ArrayList<KmerPathNode>(2);
		}
		if (to.prevList == null) {
			to.prevList = new ArrayList<KmerPathNode>(2);
		}
		from.nextList.add(to);
		to.prevList.add(from);
		if (from.nextList.size() > 1 && NEXT_SORT_ORDER.compare(from.nextList.get(from.nextList.size() - 2), from.nextList.get(from.nextList.size() - 1)) > 0) {
			from.edgesSorted = false;
		}
		if (to.prevList.size() > 1 && PREV_SORT_ORDER.compare(to.prevList.get(to.prevList.size() - 2), to.prevList.get(to.prevList.size() - 1)) > 0) {
			to.edgesSorted = false;
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			from.sanityCheck();
			to.sanityCheck();
		}
	}
	/**
	 * Replaces all edges coming to/from the source node, with edges to/from the target node
	 * @param toRemoveFrom source node to remove all edges from
	 * @param toAddTo target node to replace with
	 */
	private static void replaceEdges(KmerPathNode toRemoveFrom, KmerPathNode toAddTo) {
		if (toRemoveFrom.nextList != null) {
			if (toAddTo.nextList == null) {
				toAddTo.nextList = new ArrayList<KmerPathNode>(toRemoveFrom.nextList.size() + 2);
			}
			for (KmerPathNode n : toRemoveFrom.nextList) {
				assert(n.prevList != null);
				replaceUnique(n.prevList, toRemoveFrom, toAddTo);
				if (!toAddTo.nextList.contains(n)) {
					toAddTo.nextList.add(n);
				}
				n.edgesSorted = false;
			}
			toRemoveFrom.nextList = null;
			toAddTo.edgesSorted = false;
		}
		if (toRemoveFrom.prevList != null) {
			if (toAddTo.prevList == null) {
				toAddTo.prevList = new ArrayList<KmerPathNode>(toRemoveFrom.prevList.size() + 2);
			}
			for (KmerPathNode n : toRemoveFrom.prevList) {
				assert(n.nextList != null);
				replaceUnique(n.nextList, toRemoveFrom, toAddTo);
				if (!toAddTo.prevList.contains(n)) {
					toAddTo.prevList.add(n);
				}
				n.edgesSorted = false;
			}
			toRemoveFrom.prevList = null;
			toAddTo.edgesSorted = false;
		}
	}
	private void ensureEdgesSorted() {
		if (!edgesSorted) {
			if (nextList != null) {
				Collections.sort(nextList, NEXT_SORT_ORDER);
			}
			if (prevList != null) {
				Collections.sort(prevList, PREV_SORT_ORDER);
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
				start,
				end,
				reference,
				this.weight);
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
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			this.sanityCheck();
			split.sanityCheck();
		}
		return split;
	}
	private static void replaceFirst(List<KmerPathNode> list, KmerPathNode existing, KmerPathNode replaceWith) {
		if (list == null) return;
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i) == existing) { 
				list.set(i, replaceWith);
				return;
			}
		}
		throw new IllegalStateException("Could not replace non-existant element " + existing.toString());
	}
	/**
	 * Replaces the given element with the given replacement, ensuring that the replacement element is not
	 * duplicated.
	 * @param list
	 * @param existing
	 * @param replaceWith
	 */
	private static void replaceUnique(List<KmerPathNode> list, KmerPathNode existing, KmerPathNode replaceWith) {
		if (list == null) return;
		int existingOffset = -1;
		int replaceWithOffset = -1;
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i) == existing) {
				existingOffset = i;
			} else if (list.get(i) == replaceWith) {
				replaceWithOffset = i;
			}
		}
		if (existingOffset < 0) {
			throw new IllegalStateException("Could not replace non-existant element " + existing.toString());
		}
		if (replaceWithOffset >= 0) {
			// already in list
			list.remove(existingOffset);
		} else {
			list.set(existingOffset, replaceWith);
		}
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
		KmerPathNode split = new KmerPathNode(kmers, start, newStartPosition - 1, reference, totalWeight, weight);
		this.start = newStartPosition;
		
		if (nextList != null) {
			ArrayList<KmerPathNode> newNextThis = new ArrayList<KmerPathNode>(nextList.size());
			ArrayList<KmerPathNode> newNextSplit = new ArrayList<KmerPathNode>(nextList.size());
			for (KmerPathNode adj : nextList) {
				if (IntervalUtil.overlapsClosed(this.startPosition() + 1, this.endPosition() + 1, adj.startPosition(0), adj.endPosition(0))) {
					newNextThis.add(adj);
				}
				if (IntervalUtil.overlapsClosed(split.startPosition() + 1, split.endPosition() + 1, adj.startPosition(0), adj.endPosition(0))) {
					newNextSplit.add(adj);
				}
			}
			this.nextList = newNextThis;
			split.nextList = newNextSplit;
		}
		if (prevList != null) {
			ArrayList<KmerPathNode> newPrevThis = new ArrayList<KmerPathNode>(prevList.size());
			ArrayList<KmerPathNode> newPrevSplit = new ArrayList<KmerPathNode>(prevList.size());
			for (KmerPathNode adj : prevList) {
				if (IntervalUtil.overlapsClosed(adj.startPosition() + 1, adj.endPosition() + 1, this.startPosition(0), this.endPosition(0))) {
					newPrevThis.add(adj);
				}
				if (IntervalUtil.overlapsClosed(adj.startPosition() + 1, adj.endPosition() + 1, split.startPosition(0), split.endPosition(0))) {
					newPrevSplit.add(adj);
				}
			}
			this.prevList = newPrevThis;
			split.prevList = newPrevSplit;
		}
		split.edgesSorted = this.edgesSorted; // edge list order retained
		this.versionId++;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			this.sanityCheck();
			split.sanityCheck();
		}
		return split;
	}
	public boolean sanityCheck(int k, int maxSupportWidth, int maxPathLength) {
		sanityCheck();
		assert(length() <= maxPathLength);
		assert(end - start <= maxSupportWidth);
		for (int i = 1; i < length(); i++) {
			assert(KmerEncodingHelper.isNext(k, kmers.getLong(i - 1), kmers.getLong(i)));
		}
		assert(sumWeights(weight) == totalWeight);
		if (nextList != null) {
			for (KmerPathNode next : nextList) {
				assert(KmerEncodingHelper.isNext(k, kmer(), next.kmer(0)));
			}
		}
		if (prevList != null) {
			for (KmerPathNode prev : prevList) {
				assert(KmerEncodingHelper.isNext(k, prev.kmer(), kmer(0)));
			}
		}
		return true;
	}
	public boolean sanityCheck() {
		assert(isValid());
		assert(start <= end);
		assert(totalWeight > 0);
		assert(kmers.size() == length());
		assert(weight.size() == length());
		assert(sumWeights(weight) == totalWeight);
		if (nextList != null) {
			for (KmerPathNode next : nextList) {
				assert(IntervalUtil.overlapsClosed(startPosition() + 1, endPosition() + 1, next.startPosition(0), next.endPosition(0)));
				assert(next.isValid());
				assert(next.prev().contains(this));
			}
			if (edgesSorted) {
				assert(NEXT_SORT_ORDER.isOrdered(nextList));
			}
		}
		if (prevList != null) {
			for (KmerPathNode prev : prevList) {
				assert(IntervalUtil.overlapsClosed(startPosition(0) - 1, endPosition(0) - 1, prev.startPosition(), prev.endPosition()));
				assert(prev.isValid());
				assert(prev.next().contains(this));
			}
			if (edgesSorted) {
				assert(PREV_SORT_ORDER.isOrdered(prevList));
			}
		}
		assert(EMPTY_KMER_LIST != null && EMPTY_KMER_LIST.size() == 0); // fastutil doesn't have ImmutableList wrappers
		return true;
	}
	@Override
	public String toString() {
		if (!isValid()) {
			return String.format("[%d-%d]%s %dw (INVALID) ", startPosition(0), endPosition(0), isReference() ? "R" : " ", weight());
		}
		StringBuilder sb = new StringBuilder(String.format("[%d-%d]%s %dw ", startPosition(0), endPosition(0), isReference() ? "R" : " ", weight()));
		sb.append(KmerEncodingHelper.toApproximateString(kmer(0)));
		for (int i = 1; i < 64 && i < length(); i++) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(kmer(i)));
		}
		if (length() >= 64) {
			sb.append("...");
		}
		sb.append(String.format("(%d)", length()));
		sb.append('\n');
		return sb.toString();
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
}
