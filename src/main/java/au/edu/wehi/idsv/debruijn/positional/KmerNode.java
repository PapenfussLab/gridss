package au.edu.wehi.idsv.debruijn.positional;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public abstract class KmerNode {
	public abstract long kmer();
	public abstract int startPosition();
	public abstract int endPosition();
	public abstract int weight();
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
}
