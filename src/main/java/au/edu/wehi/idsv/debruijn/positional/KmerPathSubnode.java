package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntHeapPriorityQueue;

import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

public class KmerPathSubnode implements DeBruijnSequenceGraphNode {
	private final KmerPathNode n;
	private final int start;
	private final int end;
	public KmerPathSubnode(KmerPathNode n, int start, int end) {
		assert(n.firstKmerStart() <= start);
		assert(n.firstKmerEnd() >= end);
		assert(end - start >= 0);
		this.n = n;
		this.start = start;
		this.end = end;
	}
	public KmerPathSubnode(KmerPathNode node) {
		this(node, node.firstKmerStart(), node.firstKmerEnd());
	}
	public KmerPathNode node() { return n; }
	public long firstKmer() { return n.firstKmer(); }
	public int firstKmerStartPosition() { return start; }
	public int firstKmerEndPosition() { return end; }
	public int startPosition() { return firstKmerStartPosition() + length() - 1; }
	public int endPosition() { return firstKmerEndPosition() + length() - 1; }
	public int width() { return end - start + 1; }
	public int length() { return n.length(); }
	public int weight() { return n.weight(); }
	public int weight(int offset) { return n.weight(offset); }
	public long kmer(int offset) { return n.kmer(offset); }
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
		int targetStart = start + n.length();
		int targetEnd = end + n.length();
		for (KmerPathNode pn : n.next()) {
			int pnStart = pn.firstKmerStart();
			int pnEnd = pn.firstKmerEnd();
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
		int targetStartLastKmer = start - 1;
		int targetEndLastKmer = end - 1;
		for (KmerPathNode pn : n.prev()) {
			int pnStartLastKmer = pn.lastStart();
			int pnEndLastKmer = pn.lastEnd();
			if (pnEndLastKmer < targetStartLastKmer) {
				continue;
			} else if (pnStartLastKmer > targetEndLastKmer) {
				break;
			} else {
				adj.add(new KmerPathSubnode(pn,
						Math.max(targetStartLastKmer, pnStartLastKmer) - pn.length() + 1,
						Math.min(targetEndLastKmer, pnEndLastKmer) - pn.length() + 1));
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
		return pathRangesOfDegree(degree, true);
	}
	/**
	 * Returns the positions for which this node has the given out-degree 
	 * @param degree
	 * @return
	 */
	public RangeSet<Integer> prevPathRangesOfDegree(int degree) {
		return pathRangesOfDegree(degree, false);
	}
	/**
	 * Returns the positions for which this node has the given degree 
	 * @param degree degree expected
	 * @param outEdges calculate out edges
	 * @return
	 */
	private RangeSet<Integer> pathRangesOfDegree(int degree, boolean outEdges) {
		int start = firstKmerStartPosition();
		final int scopeEnd = firstKmerEndPosition();
		IntHeapPriorityQueue unprocessedEndsAt = new IntHeapPriorityQueue();
		IntHeapPriorityQueue unprocessedStartsAt = new IntHeapPriorityQueue();
		if (outEdges) {
			for (KmerPathSubnode adj : next()) {
				KmerPathSubnode thisRestricted = givenNext(adj);
				unprocessedStartsAt.enqueue(thisRestricted.firstKmerStartPosition());
				unprocessedEndsAt.enqueue(thisRestricted.firstKmerEndPosition());
			}
		} else {
			for (KmerPathSubnode adj : prev()) {
				KmerPathSubnode thisRestricted = givenPrev(adj);
				unprocessedStartsAt.enqueue(thisRestricted.firstKmerStartPosition());
				unprocessedEndsAt.enqueue(thisRestricted.firstKmerEndPosition());
			}
		}
		RangeSet<Integer> result = TreeRangeSet.create();
		int activeCount = 0;
		while (start <= scopeEnd) {
			// fall out of scope
			while (!unprocessedEndsAt.isEmpty() && unprocessedEndsAt.firstInt() < start) {
				unprocessedEndsAt.dequeueInt();
				activeCount--;
			}
			while (!unprocessedStartsAt.isEmpty() && unprocessedStartsAt.firstInt() <= start) {
				unprocessedStartsAt.dequeueInt();
				activeCount++;
			}
			int end = scopeEnd;
			if (!unprocessedEndsAt.isEmpty()) {
				end = Math.min(end, unprocessedEndsAt.firstInt());
			}
			if (!unprocessedStartsAt.isEmpty()) {
				end = Math.min(end, unprocessedStartsAt.firstInt() - 1);
			}
			// matches the expected number of adjacent nodes
			if (activeCount == degree) {
				result.add(Range.closed(start, end));
			}
			start = end + 1;
		}
		return result;
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s", firstKmerStartPosition(), firstKmerEndPosition(), n);
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
