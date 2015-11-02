package au.edu.wehi.idsv.debruijn.positional;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

import it.unimi.dsi.fastutil.Hash;

/**
 * KmerNode utility classes
 * @author Daniel Cameron
 *
 */
public class KmerNodeUtil {
	public static final Ordering<KmerNode> ByFirstStart = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.firstStart(), right.firstStart());
		}
	};
	public static final Ordering<KmerNode> ByFirstEnd = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.firstEnd(), right.firstEnd());
		}
	};
	public static final Ordering<KmerNode> ByFirstStartKmer = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.firstStart(), right.firstStart())
					.compare(left.firstKmer(), right.firstKmer())
					.result();
		}
	};
	public static final Ordering<KmerNode> ByFirstEndKmer = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.firstEnd(), right.firstEnd())
					.compare(left.firstKmer(), right.firstKmer())
					.result();
		}
	};
	/**
	 * Ordering usable as a graph node key (since duplicate nodes should not occur)  
	 */
	public static final Ordering<KmerNode> ByFirstStartEndKmerReference = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.firstStart(), right.firstStart())
					.compare(left.firstEnd(), right.firstEnd())
					.compare(left.firstKmer(), right.firstKmer())
					.compareTrueFirst(left.isReference(), right.isReference())
					.result();
		}
	};
	public static final Ordering<KmerNode> ByLastStart = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.lastStart(), right.lastStart());
		}
	};
	public static final Ordering<KmerNode> ByLastEnd = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.lastEnd(), right.lastEnd());
		}
	};
	public static final Ordering<KmerNode> ByLastEndKmer = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.lastEnd(), right.lastEnd())
					.compare(left.lastKmer(), right.lastKmer())
					.result();
		}
	};
	public static final Ordering<KmerNode> ByLastStartEndKmerReferenceWeight = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.lastStart(), right.lastStart())
					.compare(left.lastEnd(), right.lastEnd())
					.compare(left.lastKmer(), right.lastKmer())
					.compareTrueFirst(left.isReference(), right.isReference())
					.compare(left.weight(), right.weight())
					.result();
		}
	};
	/**
	 * Ordering usable as a graph node key (since duplicate nodes should not occur)  
	 */
	public static final Ordering<KmerNode> ByLastEndStartKmerReference = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return ComparisonChain.start()
					.compare(left.lastEnd(), right.lastEnd())
					.compare(left.lastStart(), right.lastStart())
					.compare(left.lastKmer(), right.lastKmer())
					.compareTrueFirst(left.isReference(), right.isReference())
					.result();
		}
	};
	public static final Ordering<KmerNode> ByWeight = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.weight(), right.weight());
		}
	};
	public static final Ordering<KmerNode> ByWeightDesc = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(right.weight(), left.weight());
		}
	};
	public static class HashByLastEndKmer<T extends KmerNode> implements Hash.Strategy<T> {
		@Override
		public int hashCode(T node) {
			int result = (int) (node.lastKmer() ^ (node.lastKmer() >>> 32));
			result ^= node.lastEnd();
			return result;
		}
		@Override
		public boolean equals(T a, T b) {
			return a.lastEnd() == b.lastEnd()
					&& a.lastKmer() == b.lastKmer();
		}
	}
	public static class HashByLastStartKmer<T extends KmerNode> implements Hash.Strategy<T> {
		@Override
		public int hashCode(T node) {
			int result = (int) (node.lastKmer() ^ (node.lastKmer() >>> 32));
			result ^= node.lastStart();
			return result;
		}
		@Override
		public boolean equals(T a, T b) {
			return a.lastStart() == b.lastStart()
					&& a.lastKmer() == b.lastKmer();
		}
	}
}
