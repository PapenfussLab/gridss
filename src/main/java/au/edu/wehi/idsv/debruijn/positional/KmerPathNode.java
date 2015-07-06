package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.Hash;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.util.CollectionUtil;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class KmerPathNode implements KmerNode, DeBruijnSequenceGraphNode {
	private static final LongArrayList EMPTY_KMER_LIST = new LongArrayList();
	private static final IntArrayList EMPTY_OFFSET_LIST = new IntArrayList();
	private static final List<KmerPathNode> EMPTY_EDGE_LIST = ImmutableList.of();
	private static final Ordering<KmerNode> NEXT_SORT_ORDER = KmerNodeUtil.ByFirstStart;
	private static final Ordering<KmerNode> PREV_SORT_ORDER = KmerNodeUtil.ByLastStart;
	private LongArrayList kmers; // FIXME: replace with 2-bit encoding of kmer sequence
	private LongArrayList additionalKmers = null;
	private IntArrayList additionalKmerOffsets = null;
	private IntArrayList weight;
	private int totalWeight;
	private int start;
	private int end;
	private boolean reference;
	/**
	 * Edges
	 */
	private ArrayList<KmerPathNode> nextList = null;
	private ArrayList<KmerPathNode> prevList = null;
	private boolean edgesSorted = true;
	/**
	 * Final kmer in path graph
	 */
	public long lastKmer() { return kmer(length() - 1); }
	public long firstKmer() { return kmer(0); }
	/**
	 * First possible position of final kmer
	 */
	public int lastStart() { return startPosition(length() - 1); }
	/**
	 * Last possible position of final kmer
	 */
	public int lastEnd() { return endPosition(length() - 1); }
	public int firstStart() { return start; }
	public int firstEnd() { return end; }
	public long kmer(int offset) { return kmers.getLong(offset); }
	public int startPosition(int offset) { return start + offset; }
	public int endPosition(int offset) { return end + offset; }
	public int weight() { return totalWeight; }
	public LongArrayList pathKmers() { return kmers; }
	public IntArrayList pathWeights() { return weight; }
	@Override
	public int weight(int offset) {
		return weight.getInt(offset);
	}
	public boolean isReference() { return reference; }
	public int length() { return kmers.size(); }
	public int width() { return end - start + 1; }
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
	public IntArrayList collapsedKmerOffsets()
	{
		return additionalKmerOffsets != null ? additionalKmerOffsets : EMPTY_OFFSET_LIST;
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
		this(node.lastKmer(), node.lastStart(), node.lastEnd(), node.isReference(), node.weight());
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
		assert(node.lastStart() == lastStart() + 1);
		assert(node.lastEnd() == lastEnd() + 1);
		assert(node.isReference() == isReference());
		assert(nextList == null || nextList.size() == 0);
		kmers.add(node.lastKmer());
		weight.add(node.weight());
		totalWeight += node.weight();
		reference |= node.isReference();
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			sanityCheck();
		}
	}
	/**
	 * Adds the given node to the front of this one
	 * @param node
	 */
	public void prepend(KmerPathNode node) {
		assert(firstStart() == node.lastStart() + 1);
		assert(firstEnd() == node.lastEnd() + 1);
		assert(isReference() == node.isReference());
		assert(node.nextList != null);
		assert(node.nextList.size() == 1);
		assert(node.nextList.get(0) == this);
		assert(prevList != null);
		assert(prevList.size() == 1);
		assert(prevList.get(0) == node);
		int nodeLength = node.length();
		node.kmers.addAll(kmers);
		kmers = node.kmers;
		node.weight.addAll(weight);
		weight = node.weight;
		totalWeight += node.totalWeight;
		reference |= node.reference;
		if (additionalKmerOffsets != null) {
			// shift over the additional kmers offsets to their new values
			for (int i = 0; i < additionalKmerOffsets.size(); i++) {
				additionalKmerOffsets.set(i, additionalKmerOffsets.getInt(i) + nodeLength);
			}
		}
		if (node.additionalKmers != null) {
			if (additionalKmers == null) {
				additionalKmers = new LongArrayList(node.additionalKmers.size());
				additionalKmerOffsets = new IntArrayList(node.additionalKmers.size());
			}
			additionalKmers.addAll(node.additionalKmers);
			additionalKmerOffsets.addAll(node.additionalKmerOffsets);
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
				&& weight.equals(node.weight)
				&& hasSameCollapsedKmers(node);
	}
	private boolean hasSameCollapsedKmers(KmerPathNode node) {
		if (collapsedKmers().size() != node.collapsedKmers().size()) return false;
		for (int i = 0; i < node.collapsedKmers().size(); i++) {
			if (!containsCollapsedKmer(node.additionalKmers.getLong(i), node.additionalKmerOffsets.getInt(i))) {
				return false;
			}
		}
		return true;
	}
	private boolean containsCollapsedKmer(long kmer, int offset) {
		if (additionalKmers == null) return false;
		for (int i = 0; i < additionalKmers.size(); i++) {
			if (additionalKmers.getLong(i) == kmer && additionalKmerOffsets.getInt(i) == offset) {
				return true;
			}
		}
		return false;
	}
	/**
	 * Merges the node covering the adjacent interval with matching kmers
	 * immediately preceeding this node in genomic 
	 * @param node 
	 */
	public void coaleseBeforeAdjacent(KmerPathNode node) {
		assert(canCoaleseBeforeAdjacent(node));
		// start position of this about to change so all neighbours will now be unsorted
		invalidateNeighbourEdgeSorting();
		replaceEdges(node, this);
		start = node.start;
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
		if (toMerge == this) return;
		assert(toMerge.firstStart() == firstStart());
		assert(toMerge.firstEnd() == firstEnd());
		assert(toMerge.length() == length());
		reference |= toMerge.reference;
		if (additionalKmers == null) {
			additionalKmers = toMerge.additionalKmers;
			additionalKmerOffsets = toMerge.additionalKmerOffsets;
		} else if (toMerge.additionalKmers != null) {
			additionalKmers.addAll(toMerge.additionalKmers);
			additionalKmerOffsets.addAll(toMerge.additionalKmerOffsets);
		}
		if (additionalKmers == null) {
			additionalKmers = toMerge.kmers;
			additionalKmerOffsets = new IntArrayList(toMerge.kmers.size());
		} else {
			additionalKmers.addAll(toMerge.kmers);
		}
		for (int i = 0; i < toMerge.length(); i++) {
			additionalKmerOffsets.add(i);
		}
		totalWeight += toMerge.totalWeight;
		for (int i = 0; i < weight.size(); i++) {
			weight.set(i, weight.get(i) + toMerge.weight.get(i));
		}
		replaceEdges(toMerge, this);
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			sanityCheck();
		}
		toMerge.invalidate();
	}
	/**
	 * Indicates that this node has been transformed or deleted
	 *  
	 *  Edges 
	 */
	public void invalidate() {
		assert(nextList == null || nextList.size() == 0);
		assert(prevList == null || prevList.size() == 0);
		kmers = null;
		weight = null;
		nextList = null;
		prevList = null;
		additionalKmers = null;
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
		PriorityQueue<KmerPathNode> active = new PriorityQueue<KmerPathNode>(4, KmerNodeUtil.ByFirstEnd);
		ensureEdgesSorted();
		if (nextList == null) {
			subnodes.add(new KmerPathSubnode(this));
		} else {
			int activeStart = start;
			int offset = 0;
			while (activeStart <= end) {
				while (offset < nextList.size() && nextList.get(offset).firstStart() - length() <= activeStart) {
					// add successors that should now be in scope
					active.add(nextList.get(offset++));
				}
				while (!active.isEmpty() && active.peek().firstEnd() - length() < activeStart) {
					// remove out of scope successors
					active.poll();
				}
				int endPos = end;
				if (offset < nextList.size()) {
					int nextStart = nextList.get(offset).firstStart() - length();
					endPos = Math.min(endPos, nextStart - 1);
				}
				if (!active.isEmpty()) {
					int nextEnd = active.peek().firstEnd() - length();
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
		PriorityQueue<KmerPathNode> active = new PriorityQueue<KmerPathNode>(4, KmerNodeUtil.ByFirstEnd);
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
				while (!active.isEmpty() && active.peek().lastEnd() + 1 < activeStart) {
					// remove out of scope successors
					active.poll();
				}
				int endPos = end;
				if (offset < prevList.size()) {
					int nextStart = prevList.get(offset).lastStart() + 1;
					endPos = Math.min(endPos, nextStart - 1);
				}
				if (!active.isEmpty()) {
					int nextEnd = active.peek().lastEnd() + 1;
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
		addEdgeImpl(from, to);
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			from.sanityCheck();
			to.sanityCheck();
		}
	}
	private static void addEdgeImpl(KmerPathNode from, KmerPathNode to) {
		assert(from.isValid());
		assert(to.isValid());
		assert(from.nextList == null || !CollectionUtil.containsByReference(from.nextList, to));
		assert(to.prevList == null || !CollectionUtil.containsByReference(to.prevList, from));
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
	}
	/**
	 * Replaces all edges coming to/from the source node, with edges to/from the target node
	 * @param toRemoveFrom source node to remove all edges from
	 * @param toAddTo target node to replace with
	 */
	private static void replaceEdges(KmerPathNode toRemoveFrom, KmerPathNode toAddTo) {
		assert(toRemoveFrom != toAddTo);
		assert(toRemoveFrom != null && toRemoveFrom.isValid());
		assert(toAddTo != null && toAddTo.isValid());
		if (toRemoveFrom.nextList != null) {
			if (toAddTo.nextList == null) {
				toAddTo.nextList = new ArrayList<KmerPathNode>(toRemoveFrom.nextList.size() + 2);
			}
			for (KmerPathNode n : toRemoveFrom.nextList) {
				assert(n.prevList != null);
				boolean couldRemove = CollectionUtil.removeByReference(n.prevList, toRemoveFrom);
				assert(couldRemove);
				if (CollectionUtil.addUniqueByReference(n.prevList, toAddTo)) {
					n.edgesSorted = false;
				}
				if (CollectionUtil.addUniqueByReference(toAddTo.nextList, n)) {
					toAddTo.edgesSorted = false;
				}
			}
			toRemoveFrom.nextList = null;
		}
		if (toRemoveFrom.prevList != null) {
			if (toAddTo.prevList == null) {
				toAddTo.prevList = new ArrayList<KmerPathNode>(toRemoveFrom.prevList.size() + 2);
			}
			for (KmerPathNode n : toRemoveFrom.prevList) {
				assert(n.nextList != null);
				boolean couldRemove = CollectionUtil.removeByReference(n.nextList, toRemoveFrom);
				assert(couldRemove);
				if (CollectionUtil.addUniqueByReference(n.nextList, toAddTo)) {
					n.edgesSorted = false;
				}
				if (CollectionUtil.addUniqueByReference(toAddTo.prevList, n)) {
					toAddTo.edgesSorted = false;
				}
			}
			toRemoveFrom.prevList = null;
		}
		// invalidate neighbours of our target
	}
	/**
	 * Invalidate the edge sorting of all neighbours.
	 * This should be called when the node start or end interval changes
	 */
	private void invalidateNeighbourEdgeSorting() {
		if (nextList != null) {
			for (KmerPathNode n : nextList) {
				n.edgesSorted = false;
			}
		}
		if (prevList != null) {
			for (KmerPathNode n : prevList) {
				n.edgesSorted = false;
			}
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
		split.prevList = this.prevList;
		split.edgesSorted = this.edgesSorted; // edge sorting remains unchanged for all nodes
		this.prevList = null;
		if (split.prevList != null) {
			for (KmerPathNode n : split.prevList) {
				replaceFirst(n.nextList, this, split);
			}
		}
		// outgoing edges remain unchanged since our end kmer is unchanged
		this.kmers = kmerSecond;
		this.weight = weightSecond;
		this.totalWeight -= split.weight();
		this.start += firstNodeLength;
		this.end += firstNodeLength;
		addEdgeImpl(split, this);
		if (this.additionalKmers != null) {
			LongArrayList nodeKmers = new LongArrayList();
			IntArrayList nodeOffsets = new IntArrayList();
			LongArrayList splitKmers = new LongArrayList();
			IntArrayList splitOffsets = new IntArrayList();
			for (int i = 0; i < additionalKmerOffsets.size(); i++) {
				if (additionalKmerOffsets.getInt(i) < firstNodeLength) {
					splitKmers.add(additionalKmers.getLong(i));
					splitOffsets.add(additionalKmerOffsets.getInt(i));
				} else {
					nodeKmers.add(additionalKmers.getLong(i));
					nodeOffsets.add(additionalKmerOffsets.getInt(i) - firstNodeLength);
				}
			}
			this.additionalKmers = nodeKmers;
			this.additionalKmerOffsets = nodeOffsets;
			split.additionalKmers = splitKmers;
			split.additionalKmerOffsets = splitOffsets;
		}
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
	 * Splits at the given internal start position
	 * 
	 * @param internalStart start position to split at
	 * @return new node ending immediately before the given position
	 */
	public KmerPathNode splitAtStartPosition(int newStartPosition) {
		assert(newStartPosition > start);
		assert(newStartPosition <= end);
		KmerPathNode split = new KmerPathNode(kmers, start, newStartPosition - 1, reference, totalWeight, weight);
		this.start = newStartPosition;
		if (nextList != null) {
			ArrayList<KmerPathNode> newNextThis = new ArrayList<KmerPathNode>(nextList.size());
			ArrayList<KmerPathNode> newNextSplit = new ArrayList<KmerPathNode>(nextList.size());
			for (KmerPathNode adj : nextList) {
				if (IntervalUtil.overlapsClosed(this.lastStart() + 1, this.lastEnd() + 1, adj.firstStart(), adj.firstEnd())) {
					newNextThis.add(adj);
					// start position of this has changed so that might have invalidated the sort order of the neighbour
					adj.edgesSorted = false;
				} else {
					boolean removed = CollectionUtil.removeByReference(adj.prevList, this);
					assert(removed);
				}
				if (IntervalUtil.overlapsClosed(split.lastStart() + 1, split.lastEnd() + 1, adj.firstStart(), adj.firstEnd())) {
					newNextSplit.add(adj);
					adj.prevList.add(split);
					adj.edgesSorted = false;
				}
			}
			this.nextList = newNextThis;
			split.nextList = newNextSplit;
		}
		if (prevList != null) {
			ArrayList<KmerPathNode> newPrevThis = new ArrayList<KmerPathNode>(prevList.size());
			ArrayList<KmerPathNode> newPrevSplit = new ArrayList<KmerPathNode>(prevList.size());
			for (KmerPathNode adj : prevList) {
				if (IntervalUtil.overlapsClosed(adj.lastStart() + 1, adj.lastEnd() + 1, this.firstStart(), this.firstEnd())) {
					newPrevThis.add(adj);
					// start position of this has changed so that might have invalidated the sort order of the neighbour
					adj.edgesSorted = false;
				} else {
					boolean removed = CollectionUtil.removeByReference(adj.nextList, this);
					assert(removed);
				}
				if (IntervalUtil.overlapsClosed(adj.lastStart() + 1, adj.lastEnd() + 1, split.firstStart(), split.firstEnd())) {
					newPrevSplit.add(adj);
					adj.nextList.add(split);
					adj.edgesSorted = false;
				}
			}
			this.prevList = newPrevThis;
			split.prevList = newPrevSplit;
		}
		split.edgesSorted = this.edgesSorted; // edge list order retained
		if (this.additionalKmers != null) {
			split.additionalKmers = new LongArrayList(this.additionalKmers);
			split.additionalKmerOffsets = new IntArrayList(this.additionalKmerOffsets);
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			this.sanityCheck();
			split.sanityCheck();
		}
		return split;
	}
	@Override
	public String toString() {
		if (!isValid()) {
			return String.format("[%d-%d]%s %dw (INVALID) ", firstStart(), firstEnd(), isReference() ? "R" : " ", weight());
		}
		StringBuilder sb = new StringBuilder(String.format("[%d-%d]%s %dw ", firstStart(), firstEnd(), isReference() ? "R" : " ", weight()));
		sb.append(KmerEncodingHelper.toApproximateString(firstKmer()));
		sb.replace(sb.length() - 1, sb.length(), " ");
		for (int i = 0; i < length(); i++) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(kmer(i)));
		}
		sb.append(String.format(" (%d)", length()));
		sb.append('\n');
		return sb.toString();
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + ((kmers == null) ? 0 : kmers.hashCode());
		result = prime * result + (reference ? 1231 : 1237);
		result = prime * result + start;
		result = prime * result + ((weight == null) ? 0 : weight.hashCode());
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
		if (kmers == null) {
			if (other.kmers != null)
				return false;
		} else if (!kmers.equals(other.kmers))
			return false;
		if (reference != other.reference)
			return false;
		if (start != other.start)
			return false;
		if (weight == null) {
			if (other.weight != null)
				return false;
		} else if (!weight.equals(other.weight))
			return false;
		return true;
	}
	private KmerPathNode removeKmer(int offset) {
		assert(offset >= 0);
		assert(offset < length());
		if (offset > 0 && offset < length() - 1) {
			KmerPathNode split = splitAtLength(offset + 1);
			KmerPathNode additionalSplit = split.removeKmer(offset);
			assert(additionalSplit == null);
			this.prevList = null;
			split.nextList = null;
			if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
				assert(sanityCheck());
				assert(split.sanityCheck());
			}
			return split;
		}
		// start or end kmer
		if (offset == length() - 1) {
			if (nextList != null) {
				// remove final node
				for (KmerPathNode n : nextList) {
					boolean removed = CollectionUtil.removeByReference(n.prevList, this);
					assert(removed);
				}
				nextList = null;
			}
		}
		if (offset == 0) {
			if (prevList != null) {
				for (KmerPathNode n : prevList) {
					boolean removed = CollectionUtil.removeByReference(n.nextList, this);
					assert(removed);
				}
				prevList = null;
			}
			// make sure we're still starting at the same place
			// (but after we remove ourself)
			this.start++;
			this.end++;
		}
		totalWeight -= weight.getInt(offset);
		weight.remove(offset);
		kmers.remove(offset);
		if (additionalKmers != null) {
			if (length() > 0) {
				int offsetShift = offset == 0 ? 1 : 0;
				for (int i = additionalKmers.size() - 1; i >= 0; i--) {
					int currentOffset = additionalKmerOffsets.getInt(i);
					int newOffset = currentOffset - offsetShift;
					if (newOffset < 0 || newOffset >= length()) {
						additionalKmerOffsets.remove(i);
						additionalKmers.remove(i);
					} else if (newOffset != currentOffset) {
						additionalKmerOffsets.set(i, newOffset);
					}
				}
			}
		}
		if (length() == 0) {
			invalidate();
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			if (kmers != null) {
				assert(sanityCheck());
			}
		}
		return null;
	}
	/**
	 * Removes supporting evidence from the given node
	 * @param node node to remove support 
	 * @param toRemove supporting evidence to remove. Each node 
	 * @return
	 */
	public static ArrayDeque<KmerPathNode> removeWeight(KmerPathNode node, List<? extends List<? extends KmerNode>> toRemove) {
		//KmerPathNode debugClonePathNode = new KmerPathNode(node.pathKmers(), node.firstStart(), node.firstEnd(), true, node.weight);
		//ArrayList<? extends List<? extends KmerNode>> debugCloneToRemove = Lists.newArrayList(toRemove);
		int preWeight = 0;
		int deltaWeight = 0;
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			preWeight = node.weight() * node.width();
			for (int i = 0; i < toRemove.size(); i++) {
				List<? extends KmerNode> list = toRemove.get(i);
				int start = node.startPosition(i);
				int end = node.endPosition(i);
				int nodeDeltaWeight = list == null ? 0 : list.stream().mapToInt(n -> n.weight() * IntervalUtil.overlapsWidthClosed(n.lastStart(), n.lastEnd(), start, end)).sum();
				int nodeWeight = node.weight(i) * node.width();
				if (list != null) {
					for (KmerNode n : list) {
						long nkmer = n.firstKmer();
						if (nkmer != node.kmer(i)) {
							boolean foundMatchingAltKmer = false;
							for (int j = 0; j < node.collapsedKmers().size(); j++) {
								if (node.collapsedKmers().getLong(j) == nkmer) {
									if (node.collapsedKmerOffsets().getInt(j) == i) {
										// found an alt kmer that goes to the expected node
										foundMatchingAltKmer = true;
										break;
									}
								}
							}
							assert(foundMatchingAltKmer);
						}
					}
				}
				assert(nodeDeltaWeight <= nodeWeight);
				deltaWeight += nodeDeltaWeight;
			}
			assert(deltaWeight <= preWeight);
		}
		ArrayDeque<KmerPathNode> replacement = new ArrayDeque<KmerPathNode>();
		while (!toRemove.isEmpty()) {
			int index = toRemove.size() - 1;
			List<? extends KmerNode> collection = toRemove.get(index);
			toRemove.remove(index);
			if (collection != null) {
				collection.sort(KmerNodeUtil.ByLastStart);
				node = removeWeight(replacement, node, index, collection);
				if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
					assert(replacement.stream().allMatch(n -> n.isValid()));
					if (node != null) {
						assert(node.sanityCheck());
					}
				}
			}
		}
		if (node != null) {
			assert(node.isValid());
			replacement.addFirst(node);
		}
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(replacement.stream().allMatch(n -> n.isValid()));
			replacement.stream().forEach(n -> n.sanityCheck());
			int postWeight = replacement.stream().mapToInt(n -> n.weight() * n.width()).sum();
			assert(postWeight + deltaWeight == preWeight);
			assert(replacement.stream().allMatch(pn -> pn.sanityCheck()));
		}
		return replacement;
	}
	/**
	 * Removes weight from this KmerPathNode.
	 * 
	 * Weight is removed kmer-by-kmer. In the case of removal of weight from part of the
	 * defined interval, the node is split first by length, then by position, then the
	 * weight of the resultant node is reduced
	 *  
	 * @param replacementNodeList collection to add any additional KmerPathNode created
	 * @param node node to remove weight from
	 * @param offset offset of evidence to remove
	 * @param toRemove evidence to remove
	 * @return returns the node with the earliest startPosition after removing weight
	 */
	private static KmerPathNode removeWeight(ArrayDeque<KmerPathNode> outList, KmerPathNode node, int offset, List<? extends KmerNode> toRemove) {
		if (toRemove == null || toRemove.size() == 0) return node;
		PriorityQueue<KmerNode> active = new PriorityQueue<KmerNode>(4, KmerNodeUtil.ByLastEnd);
		PeekingIterator<? extends KmerNode> startIt = Iterators.peekingIterator(toRemove.iterator());
		int start = node.startPosition(offset);
		final int scopeEnd = node.endPosition(offset);
		int weightToRemove = 0;
		KmerPathNode pre = null;
		while (start <= scopeEnd) {
			// advance
			while (startIt.hasNext() && startIt.peek().lastStart() <= start) {
				KmerNode n = startIt.next();
				weightToRemove += n.weight();
				active.add(n);
			}
			while (!active.isEmpty() && active.peek().lastEnd() < start) {
				KmerNode n = active.poll();
				weightToRemove -= n.weight();
			}
			int end = scopeEnd;
			if (startIt.hasNext() && startIt.peek().lastStart() <= end) {
				end = startIt.peek().lastStart() - 1;
			}
			if (!active.isEmpty() && active.peek().lastEnd() < end) {
				end = active.peek().lastEnd();
			}
			if (start == node.startPosition(offset) && end == node.endPosition(offset)) {
				// simple case: just down-weight this kmer
				node = removeWeight(outList, node, offset, weightToRemove);
				assert(pre == null || node == null); // only one precedecessor allowed
				return pre == null ? node : pre;
			} else {
				if (node.length() != 1) {
					// split our path node up so it's 1 kmer wide
					int postOffsetLength = node.length() - offset - 1;
					int preOffsetLength = offset;
					if (postOffsetLength > 0) {
						KmerPathNode split = node;
						node = node.splitAtLength(offset + 1);
						assert(split.isValid());
						outList.addFirst(split);
					}
					if (preOffsetLength > 0) {
						pre = node.splitAtLength(preOffsetLength);
					}
					offset = 0;
				}
				// then split path so we have a single path node containing our interval of interest 
				KmerPathNode after = null;
				if (start != node.lastStart()) {
					KmerPathNode split = node.splitAtStartPosition(start); 
					assert(split.isValid());
					outList.addFirst(split);
					assert(start == node.lastStart());
				} else if (end != node.lastEnd()) {
					after = node;
					node = node.splitAtStartPosition(end + 1);
				}
				// ok, we're ready to process this node
				if (weightToRemove > 0) {
					node = removeWeight(outList, node, 0, weightToRemove);
				}
				// move on to the next interval
				if (node != null) {
					assert(node.isValid());
					outList.addFirst(node);
				}
				node = after;
			}
			start = end + 1;
			if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
				if (node != null) {
					assert(node.sanityCheck());
				}
			}
		}
		return pre;
	}
	/**
	 * Remove the given weight from the given node 
	 * @param outList 
	 * @param node
	 * @param offset
	 * @param weightToRemove
	 * @return KmerPathNode containing kmers preceeding the offset
	 */
	private static KmerPathNode removeWeight(ArrayDeque<KmerPathNode> outList, KmerPathNode node, int offset, int weightToRemove) {
		assert(node.totalWeight >= weightToRemove);
		int newWeight = node.weight.get(offset) - weightToRemove;
		if (newWeight == 0) {
			// remove entire kmer from KmerPathNode
			KmerPathNode split = node.removeKmer(offset);
			if (split != null) {
				assert(node.isValid());
				outList.addFirst(node);
				node = split;
			}
		} else {
			node.weight.set(offset, newWeight);
			node.totalWeight -= weightToRemove;
		}
		if (!node.isValid()) return null;
		if (offset == 0) {
			// node is not valid before this position
			assert(node.isValid());
			outList.addFirst(node);
			return null;
		}
		return node;
	}
	public static class HashByFirstKmerEndPositionKmer<T extends KmerPathNode> implements Hash.Strategy<T> {
		@Override
		public int hashCode(T node) {
			int result = (int) (node.lastKmer() ^ (node.lastKmer() >>> 32));
			result ^= node.firstEnd();
			return result;
		}
		@Override
		public boolean equals(T a, T b) {
			return a.firstEnd() == b.firstEnd()
					&& a.firstKmer() == b.firstKmer();
		}
	}
	public static class HashByFirstKmerStartPositionKmer<T extends KmerPathNode> implements Hash.Strategy<T> {
		@Override
		public int hashCode(T node) {
			int result = (int) (node.lastKmer() ^ (node.lastKmer() >>> 32));
			result ^= node.lastStart();
			return result;
		}
		@Override
		public boolean equals(T a, T b) {
			return a.firstStart() == b.firstStart()
					&& a.firstKmer() == b.firstKmer();
		}
	}
	public boolean sanityCheck(int k, int maxSupportWidth, int maxPathLength) {
		//sanityCheck(); // TEMPHACK
		assert(length() <= maxPathLength);
		assert(end - start <= maxSupportWidth);
		for (int i = 1; i < length(); i++) {
			assert(KmerEncodingHelper.isNext(k, kmers.getLong(i - 1), kmers.getLong(i)));
		}
		assert(sumWeights(weight) == totalWeight);
		if (nextList != null) {
			for (KmerPathNode next : nextList) {
				assert(KmerEncodingHelper.isNext(k, lastKmer(), next.firstKmer()));
			}
		}
		if (prevList != null) {
			for (KmerPathNode prev : prevList) {
				assert(KmerEncodingHelper.isNext(k, prev.lastKmer(), firstKmer()));
			}
		}
		return true;
	}
	public boolean sanityCheck() {
		//sanityCheck(25, 1000, 1000); int temphack = 42;// TEMPHACK 
		assert(isValid());
		assert(start <= end);
		assert(totalWeight > 0);
		assert(kmers.size() == length());
		assert(weight.size() == length());
		assert(sumWeights(weight) == totalWeight);
		assert(sanityCheckEdges(this, true));
		assert(EMPTY_KMER_LIST != null && EMPTY_KMER_LIST.size() == 0); // fastutil doesn't have ImmutableList wrappers
		assert((additionalKmerOffsets == null && additionalKmers == null) || (additionalKmerOffsets != null && additionalKmers != null));
		assert((additionalKmerOffsets == null && additionalKmers == null) || additionalKmerOffsets.size() == additionalKmers.size());
		//if (additionalKmers != null) {
		//	for (int i = 0; i < additionalKmers.size(); i++) {
		//		assert(additionalKmerOffsets.getInt(i) < length());
		//		long pathKmer = kmer(additionalKmerOffsets.getInt(i));
		//		long altKmer = additionalKmers.get(i);
		//		// difference is unbounded unless only collapsing bubbles and leaves
		//		assert(KmerEncodingHelper.basesDifference(KmerEncodingHelper.MAX_K, altKmer, pathKmer) < 10);
		//	}
		//}
		//if (EvidenceTracker.TEMP_HACK_CURRENT_TRACKER != null) {
		//	EvidenceTracker.TEMP_HACK_CURRENT_TRACKER.matchesExpected(new KmerPathSubnode(this));
		//}
		return true;
	}
	private static boolean sanityCheckEdges(KmerPathNode node, boolean checkNeighbours) {
		if (node.nextList != null) {
			assert(node.nextList.stream().distinct().count() == node.nextList.size());
			for (KmerPathNode next : node.nextList) {
				assert(next != null);
				assert(next.isValid());
				assert(IntervalUtil.overlapsClosed(node.lastStart() + 1, node.lastEnd() + 1, next.firstStart(), next.firstEnd()));
				assert(CollectionUtil.containsByReference(next.prevList, node));
				if (checkNeighbours) {
					assert(sanityCheckEdges(next, false));
				}
			}
			if (node.edgesSorted) {
				assert(NEXT_SORT_ORDER.isOrdered(node.nextList));
			}
		}
		if (node.prevList != null) {
			assert(node.prevList.stream().distinct().count() == node.prevList.size());
			for (KmerPathNode prev : node.prevList) {
				assert(prev != null);
				assert(prev.isValid());
				assert(IntervalUtil.overlapsClosed(node.firstStart() - 1, node.firstEnd() - 1, prev.lastStart(), prev.lastEnd()));
				assert(CollectionUtil.containsByReference(prev.nextList, node));
				if (checkNeighbours) {
					assert(sanityCheckEdges(prev, false));
				}
			}
			if (node.edgesSorted) {
				assert(PREV_SORT_ORDER.isOrdered(node.prevList));
			}
		}
		return true;
	}
	public boolean sanityCheckReachableSubgraph() {
		Set<KmerPathNode> visited = Collections.newSetFromMap(new IdentityHashMap<KmerPathNode, Boolean>());
		Set<KmerPathNode> frontier = Collections.newSetFromMap(new IdentityHashMap<KmerPathNode, Boolean>());
		frontier.add(this);
		while (!frontier.isEmpty()) {
			KmerPathNode node = frontier.iterator().next();
			frontier.remove(node);
			if (visited.contains(node)) continue;
			visited.add(node);
			frontier.addAll(node.next());
			frontier.addAll(node.prev());
			assert(node.sanityCheck());
		}
		return true;
	}
}
