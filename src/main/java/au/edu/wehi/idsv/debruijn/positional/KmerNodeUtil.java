package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.Hash;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * KmerNode utility classes
 * @author cameron.d
 *
 */
public class KmerNodeUtil {
	public static final Ordering<KmerNode> ByFirstKmerStartPosition = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.firstKmerStartPosition(), right.firstKmerStartPosition());
		}
	};
	public static final Ordering<KmerNode> ByFirstKmerEndPosition = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.firstKmerStartPosition(), right.firstKmerStartPosition());
		}
	};
	public static final Ordering<KmerNode> ByFirstKmerStartPositionKmer = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.firstKmerStartPosition(), right.firstKmerStartPosition())
					.compare(left.kmer(), right.kmer())
					.result();
		}
	};
	/**
	 * Ordering usable as a graph node key (since duplicate nodes should not occur)  
	 */
	public static final Ordering<KmerNode> ByFirstKmerStartEndKmerReference = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.firstKmerStartPosition(), right.firstKmerStartPosition())
					.compare(left.firstKmerEndPosition(), right.firstKmerEndPosition())
					.compare(left.firstKmer(), right.firstKmer())
					.compareTrueFirst(left.isReference(), right.isReference())
					.result();
		}
	};
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
	public static class HashByEndPositionKmer<T extends KmerNode> implements Hash.Strategy<T> {
		@Override
		public int hashCode(T node) {
			int result = (int) (node.kmer() ^ (node.kmer() >>> 32));
			result ^= node.endPosition();
			return result;
		}
		@Override
		public boolean equals(T a, T b) {
			return a.endPosition() == b.endPosition()
					&& a.kmer() == b.kmer();
		}
	}
	public static class HashByStartPositionKmer<T extends KmerNode> implements Hash.Strategy<T> {
		@Override
		public int hashCode(T node) {
			int result = (int) (node.kmer() ^ (node.kmer() >>> 32));
			result ^= node.startPosition();
			return result;
		}
		@Override
		public boolean equals(T a, T b) {
			return a.startPosition() == b.startPosition()
					&& a.kmer() == b.kmer();
		}
	}
}
