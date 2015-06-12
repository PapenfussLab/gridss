package au.edu.wehi.idsv.debruijn.positional;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

public interface KmerNode {
	long kmer();
	int startPosition();
	int endPosition();
	int weight();
	boolean isReference();
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
	/**
	 * Wrapper class that considers only kmer and startPosition when determining equality. 
	 * @author cameron.d
	 *
	 */
	public class KmerStartPositionEqualityWrapper {
		private final KmerNode node;
		public KmerNode node() { return node; }
		public KmerStartPositionEqualityWrapper(KmerNode node) {
			this.node = node;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (int) (node.kmer() ^ (node.kmer() >>> 32));
			result = prime * result + node.startPosition();
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof KmerNode))
				return false;
			KmerNode other = (KmerNode) obj;
			if (node.kmer() != other.kmer())
				return false;
			if (node.startPosition() != other.startPosition())
				return false;
			return true;
		}
		
	}
	public class KmerStartPositionEqualityLookup implements KmerNode {
		private final long kmer;
		private final int start;
		public KmerStartPositionEqualityLookup(
				long kmer,
				int start) {
			this.kmer = kmer;
			this.start = start;
		}
		public long kmer() { return kmer; }
		public int startPosition() { return start; }
		public int endPosition() { return 0; }
		public int weight() { return 0; }
		public boolean isReference() { return true; }
	}
	public class KmerEndPositionEqualityWrapper {
		private final KmerNode node;
		public KmerNode node() { return node; }
		public KmerEndPositionEqualityWrapper(KmerNode node) {
			this.node = node;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (int) (node.kmer() ^ (node.kmer() >>> 32));
			result = prime * result + node.endPosition();
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof KmerNode))
				return false;
			KmerNode other = (KmerNode) obj;
			if (node.kmer() != other.kmer())
				return false;
			if (node.endPosition() != other.endPosition())
				return false;
			return true;
		}
	}
	public class KmerEndPositionEqualityLookup implements KmerNode {
		private final long kmer;
		private final int end;
		public KmerEndPositionEqualityLookup(
				long kmer,
				int end) {
			this.kmer = kmer;
			this.end = end;
		}
		public long kmer() { return kmer; }
		public int startPosition() { return 0; }
		public int endPosition() { return end; }
		public int weight() { return 0; }
		public boolean isReference() { return true; }
	}
}
