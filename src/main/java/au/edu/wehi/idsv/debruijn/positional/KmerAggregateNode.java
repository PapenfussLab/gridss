package au.edu.wehi.idsv.debruijn.positional;

import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.SortedSet;
import java.util.TreeSet;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class KmerAggregateNode extends KmerNode {
	public long kmer() { return kmer; }
	public int startPosition() { return start; }
	public int endPosition() { return end; }
	public int weight() { return weight; }
	private final long kmer;
	private final int weight;
	private final int start;
	private final int end;
	public KmerAggregateNode(long kmer, int weight, int start, int end) {
		this.kmer = kmer;
		this.weight = weight;
		this.start = start;
		this.end = end;
	}
	/**
	 * Converts a collection sorted by start position to a sorted set of aggregated intervals
	 * @param sorted node set sorted by start position
	 * @return disjoint set of aggregated intervals
	 */
	public static NavigableSet<KmerAggregateNode> aggregate(long kmer, SortedSet<? extends KmerNode> sorted) {
		TreeSet<KmerAggregateNode> result = new TreeSet<KmerAggregateNode>();
		int start = 0;
		int weight = 0;
		PriorityQueue<KmerNode> active = new PriorityQueue<KmerNode>(8, KmerNode.ByEndPosition);
		for (KmerNode s : sorted) {
			assert(s.kmer() == kmer);
			assert(s.startPosition() >= start);
			while (active.peek().endPosition() < s.startPosition()) {
				// complete
				int end = active.peek().endPosition();
				result.add(new KmerAggregateNode(kmer, weight, start, end));
				while (active.peek().endPosition() == end) {
					weight -= active.poll().weight();
				}
				start = end + 1;
			}
			if (weight > 0) {
				result.add(new KmerAggregateNode(kmer, weight, start, s.startPosition() - 1));
			}
			weight += s.weight();
			active.add(s);
		}
		while (!active.isEmpty()) {
			int end = active.peek().endPosition();
			result.add(new KmerAggregateNode(kmer, weight, start, end));
			while (active.peek().endPosition() == end) {
				weight -= active.poll().weight();
			}
			start = end + 1;
		}
		assert(weight == 0);
		return result;
	}
	@Override
	public String toString() {
		return String.format("[%d-%d] w=%d, %s", startPosition(), endPosition(), weight(), KmerEncodingHelper.toApproximateString(kmer()));
	}
}
