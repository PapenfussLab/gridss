package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public interface KmerNode {
	long kmer();
	int startPosition();
	int endPosition();
	int weight();
	boolean isReference();
	default int width() { return endPosition() - startPosition() + 1; }
	default boolean overlaps(KmerNode node) {
		return kmer() == node.kmer() && IntervalUtil.overlapsClosed(startPosition(), endPosition(), node.startPosition(), node.endPosition());
	}
	public static final Ordering<KmerNode> ByStartPosition = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.startPosition(), right.startPosition());
		}
	};
	public static final Ordering<KmerNode> ByEndPosition = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.endPosition(), right.endPosition());
		}
	};
	public static final Ordering<KmerNode> ByEndPositionKmer = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.endPosition(), right.endPosition())
					.compare(left.kmer(), right.kmer())
					.result();
		}
	};
	public static final Ordering<KmerNode> ByStartEndPositionKmerReferenceWeight = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.startPosition(), right.startPosition())
					.compare(left.endPosition(), right.endPosition())
					.compare(left.kmer(), right.kmer())
					.compareTrueFirst(left.isReference(), right.isReference())
					.compare(left.weight(), right.weight())
					.result();
		}
	};
	/**
	 * Ordering usable as a graph node key (since duplicate nodes should not occur)  
	 */
	public static final Ordering<KmerNode> ByEndStartKmerReference = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.endPosition(), right.endPosition())
					.compare(left.startPosition(), right.startPosition())
					.compare(left.kmer(), right.kmer())
					.compareTrueFirst(left.isReference(), right.isReference())
					.result();
		}
	};
}
