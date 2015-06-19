package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.Hash;

/**
 * KmerNode utility classes
 * @author cameron.d
 *
 */
class KmerNodeUtil {
	static class HashByEndPositionKmer<T extends KmerNode> implements Hash.Strategy<T> {
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
	static class HashByStartPositionKmer<T extends KmerNode> implements Hash.Strategy<T> {
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
