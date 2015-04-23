package au.edu.wehi.idsv.debruijn.positional;

import java.util.Collection;
import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.TreeSet;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

public class KmerAggregateNode extends KmerNode {
	public long getKmer() { return kmer; }
	public int getStartPosition() { return start; }
	public int getEndPosition() { return end; }
	public int getWeight() { return weight; }
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
	public static NavigableSet<KmerAggregateNode> aggregate(long kmer, Collection<KmerSupportNode> sorted) {
		TreeSet<KmerAggregateNode> result = new TreeSet<KmerAggregateNode>();
		int start = 0;
		int weight = 0;
		PriorityQueue<KmerSupportNode> active = new PriorityQueue<KmerSupportNode>(8, KmerNode.ByEndPosition);
		for (KmerSupportNode s : sorted) {
			assert(s.getKmer() == kmer);
			assert(s.getStartPosition() >= start);
			while (active.peek().getEndPosition() < s.getStartPosition()) {
				// complete
				int end = active.peek().getEndPosition();
				result.add(new KmerAggregateNode(kmer, weight, start, end));
				while (active.peek().getEndPosition() == end) {
					weight -= active.poll().getWeight();
				}
				start = end + 1;
			}
			if (weight > 0) {
				result.add(new KmerAggregateNode(kmer, weight, start, s.getStartPosition() - 1));
			}
			weight += s.getWeight();
			active.add(s);
		}
		while (!active.isEmpty()) {
			int end = active.peek().getEndPosition();
			result.add(new KmerAggregateNode(kmer, weight, start, end));
			while (active.peek().getEndPosition() == end) {
				weight -= active.poll().getWeight();
			}
			start = end + 1;
		}
		assert(weight == 0);
		return result;
	}
	@Override
	public String toString() {
		return String.format("[%d-%d] w=%d, %s", getStartPosition(), getEndPosition(), getWeight(), KmerEncodingHelper.toApproximateString(getKmer()));
	}
}
