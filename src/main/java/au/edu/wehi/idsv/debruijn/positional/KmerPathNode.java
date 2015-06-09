package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * Total support for the given kmer over the given interval
 * @author cameron.d
 *
 */
public class KmerPathNode extends KmerNode {
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
	public int weight() { return weight; }
	public boolean isReference() { return reference; }
	public int length() { return kmers.size(); }
	private final List<Long> kmers;
	private int weight;
	private final int start;
	private final int end;
	private final boolean reference;
	public KmerPathNode(long kmer, int weight, int start, int end, boolean reference) {
		this.kmers = Lists.newArrayList(kmer);
		this.weight = weight;
		this.start = start;
		this.end = end;
		this.reference = reference;
	}
	public KmerPathNode(KmerNode node) {
		this(node.kmer(), node.weight(), node.startPosition(), node.endPosition(), node.isReference());
	}
	public void append(KmerAggregateNode node) {
		assert(node.startPosition() == startPosition(length() - 1) + 1);
		assert(node.endPosition() == endPosition(length() - 1) + 1);
		assert(node.isReference() == isReference());
		kmers.add(node.kmer());
		weight += node.weight();
	}
	/**
	 * Edges
	 */
	List<KmerPathNode> nextList = new ArrayList<KmerPathNode>(3);
	List<KmerPathNode> prevList = new ArrayList<KmerPathNode>(3);
	public List<KmerPathNode> next() {
		return nextList;
	}
	public List<KmerPathNode> prev() {
		return prevList;
	}
	public static void addEdge(KmerPathNode from, KmerPathNode to) {
		from.nextList.add(to);
		to.prevList.add(from);
	}
	@Override
	public String toString() {
		return String.format("[%d-%d]%s%d, %d kmers", isReference() ? "R" : " ",  startPosition(), endPosition(), weight(), length());
	}
	public static final Ordering<KmerPathNode> ByFirstKmerStartPosition = new Ordering<KmerPathNode>() {
		@Override
		public int compare(KmerPathNode left, KmerPathNode right) {
			return Ints.compare(left.startPosition(0), right.startPosition(0));
		}
	};
	public static final Ordering<KmerPathNode> ByLastKmerEndPosition = new Ordering<KmerPathNode>() {
		@Override
		public int compare(KmerPathNode left, KmerPathNode right) {
			return Ints.compare(left.startPosition(), right.startPosition());
		}
	};
}
