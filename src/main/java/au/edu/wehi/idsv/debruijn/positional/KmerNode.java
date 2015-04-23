package au.edu.wehi.idsv.debruijn.positional;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public abstract class KmerNode {
	public abstract int getStartPosition();
	public abstract int getEndPosition();
	public abstract int getWeight();
	public static final Ordering<KmerNode> ByStartPosition = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.getStartPosition(), right.getStartPosition());
		}
	};
	public static final Ordering<KmerNode> ByEndPosition = new Ordering<KmerNode>() {
		@Override
		public int compare(KmerNode left, KmerNode right) {
			return Ints.compare(left.getEndPosition(), right.getEndPosition());
		}
	};
}
