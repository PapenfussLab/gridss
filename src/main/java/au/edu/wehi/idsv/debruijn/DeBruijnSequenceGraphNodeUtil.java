package au.edu.wehi.idsv.debruijn;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.PrimitiveIterator.OfLong;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
import java.util.stream.StreamSupport;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.graph.WeightedSequenceGraphNodeUtil;


public class DeBruijnSequenceGraphNodeUtil {
	/**
	 * Returns the number of bases difference between the two paths
	 * Bases are only compared up to the length of the shortest of the
	 * paths
	 * @param pathA
	 * @param pathB
	 * @return number of bases difference between the two paths
	 */
	public static int basesDifferent(int k, Iterable<? extends DeBruijnSequenceGraphNode> pathA, Iterable<? extends DeBruijnSequenceGraphNode> pathB) {
		return basesDifferent(k, pathA, pathB, 0, 0);
	}
	/**
	 * Returns the number of bases difference between the two paths when
	 * traversing the paths backward.
	 * Bases are only compared up to the length of the shortest of the
	 * paths
	 * @param pathA
	 * @param pathB
	 * @return number of bases difference between the two paths
	 */
	public static int reverseBasesDifferent(int k, Iterable<? extends DeBruijnSequenceGraphNode> pathA, Iterable<? extends DeBruijnSequenceGraphNode> pathB) {
		int lengthA = WeightedSequenceGraphNodeUtil.nodeLength(pathA);
		int lengthB = WeightedSequenceGraphNodeUtil.nodeLength(pathB);
		int skipCountA = Math.max(0, lengthA - lengthB);
		int skipCountB = Math.max(0, lengthB - lengthA);
		// skip initial bases of the longer path
		int diff = basesDifferent(k, pathA, pathB, skipCountA, skipCountB); 
		return diff;
	}
	/**
	 * Returns the number of bases difference between the two paths
	 * @param k
	 * @param pathA
	 * @param pathB
	 * @param forwardKmerTraversal sequence to successor next kmers 
	 * @return number of bases difference
	 */
	public static int basesDifferent(int k, LongStream pathA, LongStream pathB, boolean forwardKmerTraversal) {
		OfLong itA = pathA.iterator();
		OfLong itB = pathB.iterator();
		if (!itA.hasNext() || !itB.hasNext()) return 0;
		int diff = KmerEncodingHelper.basesDifference(k, itA.nextLong(), itB.nextLong());
		while (itA.hasNext() && itB.hasNext()) {
			if (forwardKmerTraversal) {
				if (!KmerEncodingHelper.lastBaseMatches(k, itA.nextLong(), itB.nextLong())) {
					diff++;
				}
			} else {
				if (!KmerEncodingHelper.firstBaseMatches(k, itA.nextLong(), itB.nextLong())) {
					diff++;
				}
			}
		}
		return diff;
	}
	private static int basesDifferent(int k, Iterable<? extends DeBruijnSequenceGraphNode> pathA, Iterable<? extends DeBruijnSequenceGraphNode> pathB, final int initialSkipCountA, final int initialSkipCountB) {
		int skipCountA = initialSkipCountA;
		int skipCountB = initialSkipCountB;
		int diff = 0;
		Iterator<? extends DeBruijnSequenceGraphNode> itA = pathA.iterator();
		Iterator<? extends DeBruijnSequenceGraphNode> itB = pathB.iterator();
		DeBruijnSequenceGraphNode currentA = itA.hasNext() ? itA.next() : null;
		DeBruijnSequenceGraphNode currentB = itB.hasNext() ? itB.next() : null;
		int offsetA = 0;
		int offsetB = 0;
		boolean isFirstKmer = true;
		while (currentA != null && currentB != null) {
			// advance to next kmer
			while (currentA != null && offsetA >= currentA.length()) {
				offsetA = 0;
				currentA = null;
				if (itA.hasNext()) {
					currentA = itA.next();
				}
			}
			while (currentB != null && offsetB >= currentB.length()) {
				offsetB = 0;
				currentB = null;
				if (itB.hasNext()) {
					currentB = itB.next();
				}
			}
			if (currentA != null && currentB != null) {
				if (skipCountA > 0) {
					offsetA++;
					skipCountA--;
				} else if (skipCountB > 0) {
					offsetB++;
					skipCountB--;
				} else {
					// compare bases
					if (isFirstKmer) {
						diff = KmerEncodingHelper.basesDifference(k, currentA.kmer(offsetA), currentB.kmer(offsetB));
						isFirstKmer = false;
					} else if (!KmerEncodingHelper.lastBaseMatches(k, currentA.kmer(offsetA), currentB.kmer(offsetB))) {
						diff++;
					}
					offsetA++;
					offsetB++;
				}
			}
		}
		if (Defaults.SANITY_CHECK_DE_BRUIJN) {
			int streamDiff = basesDifferent(k,
					StreamSupport.stream(pathA.spliterator(), false).flatMapToLong(n -> IntStream.range(0, n.length()).mapToLong(i -> n.kmer(i))).skip(initialSkipCountA),
					StreamSupport.stream(pathB.spliterator(), false).flatMapToLong(n -> IntStream.range(0, n.length()).mapToLong(i -> n.kmer(i))).skip(initialSkipCountB),
					true);
			assert(diff == streamDiff);
		}
		return diff;
	}
	public static List<Long> asKmers(final Iterable<? extends DeBruijnSequenceGraphNode> path) {
		List<Long> kmers = new ArrayList<Long>();
		for (DeBruijnSequenceGraphNode n : path) {
			for (int i = 0; i < n.length(); i++) {
				kmers.add(n.kmer(i));
			}
		}
		return kmers;
	}
}
