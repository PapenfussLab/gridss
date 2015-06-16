package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntHeapPriorityQueue;

import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

public class KmerPathSubnode implements DeBruijnSequenceGraphNode {
	private final KmerPathNode n;
	private final int start;
	private final int end;
	public KmerPathSubnode(KmerPathNode n, int start, int end) {
		assert(n.startPosition(0) <= start);
		assert(n.endPosition(0) >= end);
		assert(end - start >= 0);
		this.n = n;
		this.start = start;
		this.end = end;
	}
	public KmerPathSubnode(KmerPathNode node) {
		this(node, node.startPosition(0), node.endPosition(0));
	}
	public KmerPathNode node() { return n; }
	public int firstKmerStartPosition() { return start; }
	public int firstKmerEndPosition() { return end; }
	public int width() { return start - end + 1; }
	@Override
	public int length() {
		return n.length();
	}
	@Override
	public int weight() {
		return n.weight();
	}
	@Override
	public int weight(int offset) {
		return n.weight(offset);
	}
	@Override
	public long kmer(int offset) {
		return n.kmer(offset);
	}
	/**
	 * Determines whether this node is a positional subset of the given node
	 * @param node node to check containment of
	 * @return true if this node is a subset of the given node, false otherwise
	 */
	public boolean isSubsetOf(KmerPathSubnode node) {
		return node.n == n
				&& node.start <= start
				&& node.end >= end;
	}
	/**
	 * Returns the subset of valid position of this node for the given next traversal 
	 * @param node next node traversed
	 * @return subset of this that is valid for the given traversal
	 */
	public KmerPathSubnode givenNext(KmerPathSubnode node) {
		int subsetStart = Math.max(firstKmerStartPosition(), node.firstKmerStartPosition() - this.length());
		int subsetEnd = Math.min(firstKmerEndPosition(), node.firstKmerEndPosition() - this.length());
		if (subsetStart == firstKmerStartPosition() && subsetEnd == firstKmerEndPosition()) {
			return this;
		} else {
			return new KmerPathSubnode(node(), subsetStart, subsetEnd);
		}
	}
	/**
	 * Returns the subset of valid position of this node for the given prev traversal 
	 * @param node prev node traversed
	 * @return subset of this that is valid for the given traversal
	 */
	public KmerPathSubnode givenPrev(KmerPathSubnode node) {
		int subsetStart = Math.max(firstKmerStartPosition(), node.firstKmerStartPosition() + node.length());
		int subsetEnd = Math.min(firstKmerEndPosition(), node.firstKmerEndPosition() + node.length());
		if (subsetStart == firstKmerStartPosition() && subsetEnd == firstKmerEndPosition()) {
			return this;
		} else {
			return new KmerPathSubnode(node(), subsetStart, subsetEnd);
		}
	}
	public List<KmerPathSubnode> next() {
		List<KmerPathSubnode> adj = new ArrayList<KmerPathSubnode>(n.next().size());
		int targetStart = start + n.length() + 1;
		int targetEnd = end + n.length() + 1;
		for (KmerPathNode pn : n.next()) {
			int pnStart = pn.startPosition(0);
			int pnEnd = pn.endPosition(0);
			// since next() is sorted, we only need to process the neighbours overlapping our interval
			if (pnEnd < targetStart) {
				continue;
			} else if (pnStart > targetEnd) {
				break;
			} else {
				adj.add(new KmerPathSubnode(pn,
						Math.max(targetStart, pnStart),
						Math.min(targetEnd, pnEnd)));
			}
		}
		return adj;
	}
	public List<KmerPathSubnode> prev() {
		List<KmerPathSubnode> adj = new ArrayList<KmerPathSubnode>(n.prev().size());
		int targetStart = start - 1;
		int targetEnd = end - 1;
		for (KmerPathNode pn : n.prev()) {
			int pnStart = pn.startPosition();
			int pnEnd = pn.endPosition();
			if (pnEnd < targetStart) {
				continue;
			} else if (pnStart > targetEnd) {
				break;
			} else {
				adj.add(new KmerPathSubnode(pn,
						Math.max(targetStart, pnStart) - pn.length(),
						Math.min(targetEnd, pnEnd) - pn.length()));
			}
		}
		return adj;
	}
	/**
	 * Returns the positions for which this node has the given out-degree 
	 * @param degree
	 * @return
	 */
	public RangeSet<Integer> nextPathRangesOfDegree(int degree) {
		return pathRangesOfDegree(degree, Iterators.peekingIterator(Iterators.transform(next().iterator(), new Function<KmerPathSubnode, KmerPathSubnode>(){
					@Override
					public KmerPathSubnode apply(KmerPathSubnode input) {
						return givenNext(input);
					}})));
	}
	/**
	 * Returns the positions for which this node has the given out-degree 
	 * @param degree
	 * @return
	 */
	public RangeSet<Integer> prevPathRangesOfDegree(int degree) {
		return pathRangesOfDegree(degree, Iterators.peekingIterator(Iterators.transform(prev().iterator(), new Function<KmerPathSubnode, KmerPathSubnode>(){
					@Override
					public KmerPathSubnode apply(KmerPathSubnode input) {
						return givenPrev(input);
					}})));
	}
	/**
	 * Returns the positions for which this node has the given degree 
	 * @param degree
	 * @return
	 */
	private RangeSet<Integer> pathRangesOfDegree(int degree, PeekingIterator<KmerPathSubnode> it) {
		RangeSet<Integer> result = TreeRangeSet.create();
		IntHeapPriorityQueue endsAt = new IntHeapPriorityQueue();
		int start = firstKmerStartPosition();
		final int scopeEnd = firstKmerEndPosition();
		while (start <= scopeEnd) {
			while (it.hasNext() && it.peek().firstKmerStartPosition() <= start) {
				endsAt.enqueue(it.peek().firstKmerStartPosition());
			}
			while (!endsAt.isEmpty() && endsAt.firstInt() < start) {
				endsAt.dequeueInt();
			}
			int end = scopeEnd;
			if (it.hasNext()) {
				end = Math.min(end, it.peek().firstKmerEndPosition() - 1);
			}
			if (!endsAt.isEmpty()) {
				end = Math.min(end, endsAt.firstInt());
			}
			// matches the expected number of adjacent nodes
			if (endsAt.size() == degree) {
				result.add(Range.closed(start, end));
			}
			start = end + 1;
		}
		return result;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + ((n == null) ? 0 : n.hashCode());
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
		KmerPathSubnode other = (KmerPathSubnode) obj;
		if (end != other.end)
			return false;
		if (n == null) {
			if (other.n != null)
				return false;
		} else if (!n.equals(other.n))
			return false;
		if (start != other.start)
			return false;
		return true;
	}
}
