package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.debruijn.DeBruijnSequenceGraphNode;
import au.edu.wehi.idsv.util.RangeUtil;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.ints.IntHeapPriorityQueue;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntPredicate;

public class KmerPathSubnode implements DeBruijnSequenceGraphNode, KmerNode {
	private static final Log log = Log.getInstance(KmerPathSubnode.class);
	public static final IntPredicate NO_EDGES = n -> n == 0;
	public static final IntPredicate SINGLE_EDGE = n -> n == 1;
	public static final IntPredicate MULTIPLE_EDGES = n -> n > 1;
	public static final IntPredicate NOT_MULTIPLE_EDGES = n -> n <= 1;
	private final KmerPathNode n;
	private final int start;
	private final int end;
	public KmerPathSubnode(KmerPathNode n, int start, int end) {
		assert(n != null);
		assert(n.firstStart() <= start);
		assert(n.firstEnd() >= end);
		assert(end - start >= 0);
		this.n = n;
		this.start = start;
		this.end = end;
	}
	public KmerPathSubnode(KmerPathNode node) {
		this(node, node.firstStart(), node.firstEnd());
	}
	public KmerPathNode node() { return n; }
	public long firstKmer() { return n.firstKmer(); }
	public long lastKmer() { return n.lastKmer(); }
	public int firstStart() { return start; }
	public int firstEnd() { return end; }
	public int lastStart() { return firstStart() + length() - 1; }
	public int lastEnd() { return firstEnd() + length() - 1; }
	public int width() { return end - start + 1; }
	public int length() { return n.length(); }
	public int weight() { return n.weight(); }
	public int weight(int offset) { return n.weight(offset); }
	public long kmer(int offset) { return n.kmer(offset); }
	public boolean isReference() { return n.isReference(); }
	/**
	 * Returns the subset of valid position of this node for the given next traversal 
	 * @param node next node traversed
	 * @return subset of this that is valid for the given traversal
	 */
	public KmerPathSubnode givenNext(KmerPathSubnode node) {
		int subsetStart = Math.max(firstStart(), node.firstStart() - this.length());
		int subsetEnd = Math.min(firstEnd(), node.firstEnd() - this.length());
		if (subsetStart == firstStart() && subsetEnd == firstEnd()) {
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
		int subsetStart = Math.max(firstStart(), node.firstStart() + node.length());
		int subsetEnd = Math.min(firstEnd(), node.firstEnd() + node.length());
		if (subsetStart == firstStart() && subsetEnd == firstEnd()) {
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
			int pnStart = pn.firstStart();
			int pnEnd = pn.firstEnd();
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
	public RangeSet<Integer> nextPathRangesOfDegree(IntPredicate expectedDegree) {
		return pathRangesOfDegree(expectedDegree, true);
	}
	/**
	 * Returns the positions for which this node has the given out-degree 
	 * @param degree
	 * @return
	 */
	public RangeSet<Integer> prevPathRangesOfDegree(IntPredicate expectedDegree) {
		return pathRangesOfDegree(expectedDegree, false);
	}
	/**
	 * Returns the position subinterval of this node for which there
	 * exist exactly the given number of adjacent nodes 
	 * @param prevDegree number of predecessor nodes
	 * @param nextDegree number of successor nodes
	 * @return Nodes with the given predecessor and successor counts
	 */
	public List<KmerPathSubnode> subnodesOfDegree(IntPredicate prevDegree, IntPredicate nextDegree) {
		RangeSet<Integer> rsNext = nextPathRangesOfDegree(nextDegree);
		if (rsNext.isEmpty()) return ImmutableList.of();
		RangeSet<Integer> rsPrev = prevPathRangesOfDegree(prevDegree);
		if (rsPrev.isEmpty()) return ImmutableList.of();
		RangeSet<Integer> rs = RangeUtil.intersect(rsNext, rsPrev);
		if (rs.isEmpty()) return ImmutableList.of();
		if (rs.encloses(Range.closed(start, end))) return ImmutableList.of(this);
		List<KmerPathSubnode> list = new ArrayList<KmerPathSubnode>(rs.asRanges().size());
		for (Range<Integer> range : rs.asRanges()) {
			list.add(new KmerPathSubnode(n, range.lowerEndpoint(), range.upperEndpoint()));
		}
		return list;
	}
	/**
	 * Returns the positions for which this node has the given degree 
	 * @param degree degree expected
	 * @param outEdges calculate out edges
	 * @return
	 */
	private RangeSet<Integer> pathRangesOfDegree(IntPredicate degreeMatchingFunction, boolean outEdges) {
		int start = firstStart();
		final int scopeEnd = firstEnd();
		IntHeapPriorityQueue unprocessedEndsAt = new IntHeapPriorityQueue();
		IntHeapPriorityQueue unprocessedStartsAt = new IntHeapPriorityQueue();
		if (outEdges) {
			for (KmerPathSubnode adj : next()) {
				KmerPathSubnode thisRestricted = givenNext(adj);
				unprocessedStartsAt.enqueue(thisRestricted.firstStart());
				unprocessedEndsAt.enqueue(thisRestricted.firstEnd());
			}
		} else {
			for (KmerPathSubnode adj : prev()) {
				KmerPathSubnode thisRestricted = givenPrev(adj);
				unprocessedStartsAt.enqueue(thisRestricted.firstStart());
				unprocessedEndsAt.enqueue(thisRestricted.firstEnd());
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
			if (degreeMatchingFunction.test(activeCount)) {
				result.add(Range.closed(start, end));
			}
			start = end + 1;
		}
		return result;
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s", firstStart(), firstEnd(), n);
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
	public boolean sanityCheck() {
		assert(n != null);
		if (!n.isValid()) {
			log.error(String.format("[%d,%d] subnode over invalid node", firstStart(), firstEnd()));
		}
		assert(n.isValid());
		assert(start >= n.firstStart());
		assert(end <= n.firstEnd());
		return true;
	}
}
