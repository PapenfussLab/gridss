package au.edu.wehi.idsv.debruijn;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
		return basesDifferent(k, pathA, pathB, skipCountA, skipCountB);
	}
	private static int basesDifferent(int k, Iterable<? extends DeBruijnSequenceGraphNode> pathA, Iterable<? extends DeBruijnSequenceGraphNode> pathB, int skipCountA, int skipCountB) {
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
