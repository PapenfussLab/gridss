package au.edu.wehi.socrates.debruijn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import org.apache.commons.lang3.ArrayUtils;

import au.edu.wehi.socrates.BreakendDirection;
import au.edu.wehi.socrates.sam.AnomolousReadAssembly;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;

public class DeBruijnReadGraph {
	private static final int MAX_QUAL_SCORE = 128 - 66;
	private final Map<Long, DeBruijnNode> kmers = Maps.newHashMap();
	private final Multimap<Long, Integer> startkmers = HashMultimap.<Long, Integer>create();
	private final int k;
	private final BreakendDirection direction;
	public DeBruijnReadGraph(int k, BreakendDirection direction) {
		this.k = k;
		this.direction = direction;
	}
	private class AnchorOffset {
		private AnchorOffset(int readOffset) {
			this.readOffset = readOffset;
		}
		private final int readOffset;
		/**
		 * kmer offset of the anchor in the read 
		 */
		public int getReadOffset() {
			return Math.max(0, readOffset);
		}
		/**
		 * Offset of breakpoint from end of the kmer
		 */
		public int getBreakpointOffset() {
			// soft-clip size is greater than our kmer length
			// since the breakpoint does not occur after a kmer, we record an offset to work out the actual breakpoint
			return Math.min(0, readOffset);
		}
	}
	private AnchorOffset getAnchoredBreakpointStartKmer(SAMRecord record, boolean anchored) {
		if (!anchored) return null;
		int clippedBases = direction == BreakendDirection.Forward ?
				record.getUnclippedEnd() - record.getAlignmentEnd() : 
				record.getAlignmentStart() - record.getUnclippedStart();
		int anchoredBreakpointStartKmer = record.getReadLength() - k - clippedBases;
		if (anchoredBreakpointStartKmer >= record.getReadLength() - k + 1) {
			throw new RuntimeException(String.format("Sanity check failure: breakpoint at %dth kmer read of length %d with k=%d", anchoredBreakpointStartKmer, record.getReadLength(), k));
		}
		return new AnchorOffset(anchoredBreakpointStartKmer);
	}
	private boolean shouldReverse(SAMRecord record, boolean anchored) {
		boolean reverseDirection = direction == BreakendDirection.Backward;
		return reverseDirection ^ reverseComplimentRequiredForPositiveStrandReadout(record, anchored);
	}
	private boolean reverseComplimentRequiredForPositiveStrandReadout(SAMRecord record, boolean anchored) {
		if (anchored && record.getReadUnmappedFlag()) throw new RuntimeException("Sanity check failure: anchored read is unmapped");
		if (anchored) return false;
		if (!record.getReadPairedFlag()) throw new RuntimeException("Sanity check failure: Cannot determine expected orientation of unanchored unpaired read");
		if (record.getMateUnmappedFlag()) throw new RuntimeException("Sanity check failure: Cannot determine expected orientation of unanchored read with unmapped mate");
		// We expect FR pairing orientation
		// Mate anchored on negative strand = need positive strand mate
		// 	-> reverse if mate+ & read+ or mate- & read- === mate XOR read
		return record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
	}
	private boolean shouldCompliment(SAMRecord record, boolean anchored) {
		return reverseComplimentRequiredForPositiveStrandReadout(record, anchored);
	}
	public void addRead(SAMRecord record, boolean anchored) {
		if (record == null) return;
		addRead(record, record.getReadBases(), record.getBaseQualities(),
				getAnchoredBreakpointStartKmer(record, anchored),
				shouldReverse(record, anchored),
				shouldCompliment(record, anchored));
	}
	private void addRead(SAMRecord record, byte[] bases, byte[] qual, AnchorOffset anchor, boolean reverse, boolean complement) {
		bases = fixReadBases(bases, reverse, complement);
		qual = fixReadQuals(qual, reverse, complement);
		int kmersTillAnchor = anchor == null ? -1 : anchor.getReadOffset();
		for (ReadKmer kmer : new ReadKmerIterable(k, bases, qual)) {
			DeBruijnNode node = kmers.get(kmer.kmer);
			if (node == null) {
				node = new DeBruijnNode();
				kmers.put(kmer.kmer, node);
			}
			node.addRead(kmer.weight, record);
			if (kmersTillAnchor-- == 0) {
				startkmers.put(kmer.kmer, anchor.getBreakpointOffset());
			}
		}
	}
	private byte[] fixReadBases(byte[] bases, boolean reverse, boolean complement) {
		if (reverse || complement) {
			bases = ArrayUtils.clone(bases);
		}
		if (reverse) {
			ArrayUtils.reverse(bases);
		}
		if (complement) {
			for (int i = 0; i < bases.length; i++) {
				bases[i] = SequenceUtil.complement(bases[i]); 
			}
		}
		return bases;
	}
	private byte[] fixReadQuals(byte[] qual, boolean reverse, boolean complement) {
		if (reverse) {
			qual = ArrayUtils.clone(qual);
			ArrayUtils.reverse(qual);
		}
		return qual;
	}
	private void removeRead(SAMRecord record, byte[] bases, byte[] qual, AnchorOffset anchor, boolean reverse, boolean complement) {
		bases = fixReadBases(bases, reverse, complement);
		qual = fixReadQuals(qual, reverse, complement);
		int kmersTillAnchor = anchor == null ? -1 : anchor.getReadOffset();
		for (ReadKmer kmer : new ReadKmerIterable(k, bases, qual)) {
			DeBruijnNode node = kmers.get(kmer.kmer);
			if (node == null) {
				throw new RuntimeException("Sanity check failure: attempt to remove kmer not in de Bruijn graph");
			}
			if (node.removeRead(kmer.weight, record)) {
				kmers.remove(kmer.kmer);
			}
			if (kmersTillAnchor-- == 0) {
				startkmers.remove(kmer.kmer, anchor.getBreakpointOffset());
			}
		}
	}
	public void removeRead(SAMRecord record, boolean anchored) {
		if (record == null) return;
		removeRead(record, record.getReadBases(), record.getBaseQualities(),
				getAnchoredBreakpointStartKmer(record, anchored),
				shouldReverse(record, anchored),
				shouldCompliment(record, anchored));
	}
	public AnomolousReadAssembly assembleVariant() {
		// debugPrint();
		return greedyAnchoredTraverse();
	}
	/**
	 * Simple greedy traversal starting from the highest weighted starting node
	 * with no repeated nodes
	 * @return
	 */
	private AnomolousReadAssembly greedyAnchoredTraverse() {
		//debugPrint();
		Long start = bestStartingPosition();
		if (start == null) return null;
		LinkedList<Long> path = greedyTraverse(start);
		AnomolousReadAssembly result = pathToAnomolousReadAssembly(path, start);
		return result;
	}
	private Long bestStartingPosition() {
		Long best = null;
		long bestScore = -1;
		for (Long state : startkmers.keySet()) {
			long weight = kmers.get(state).getWeight();
			if (weight > bestScore) {
				bestScore = weight;
				best = state;
			}
		}
		return best;
	}
	private AnomolousReadAssembly pathToAnomolousReadAssembly(LinkedList<Long> path, Long breakpointAnchor) {
		if (path == null || path.size() == 0) throw new IllegalArgumentException("Invalid path");
		int assemblyLength = path.size() + k - 1;
		int readBaseCount = 0;
		byte[] bases = KmerEncodingHelper.encodedToPicardBases(path.get(0), k);
		Set<SAMRecord> lastNodeSupport = Sets.newHashSet();
		bases = Arrays.copyOf(bases, assemblyLength);
		int offset = k - 1;
		List<Long> qual = new ArrayList<Long>(path.size());
		int softclipSize = 0;
		for (Long node : path) {
			bases[offset] = KmerEncodingHelper.lastBaseEncodedToPicardBase(node, k);
			// subtract # reads to adjust for the +1 qual introduced by ReadKmerIterable
			// to ensure positive node weights
			qual.add(this.kmers.get(node).getWeight() - this.kmers.get(node).getSupportingReads().size());
			readBaseCount += this.kmers.get(node).getSupportingReads().size();
			offset++;
			if (node == breakpointAnchor) {
				softclipSize = assemblyLength - offset;
			}
			for (SAMRecord read : this.kmers.get(node).getSupportingReads()) {
				if (lastNodeSupport.contains(read)) {
					readBaseCount++;
				} else {
					readBaseCount += k;
				}
			}
			lastNodeSupport = this.kmers.get(node).getSupportingReads();
		}
		if (!startkmers.containsEntry(breakpointAnchor, 0)) {
			// the breakpoint anchor is actually within the kmer, not at the end
			// just grab the the first offset.
			softclipSize -= startkmers.get(breakpointAnchor).iterator().next();
		}
		// pad out qualities to match the path length
		for (int i = 0; i < k - 1; i++) qual.add(qual.get(qual.size() - 1));
		byte[] quals = rescaleBaseQualities(qual);
		if (direction == BreakendDirection.Backward) {
			ArrayUtils.reverse(bases);
			ArrayUtils.reverse(quals);
		}
		return new AnomolousReadAssembly("idsvDeBruijn", bases, quals, assemblyLength - softclipSize, direction, getPathSupportingReads(path).size(), readBaseCount);
	}
	private Set<SAMRecord> getPathSupportingReads(LinkedList<Long> path) {
		Set<SAMRecord> reads = Sets.newHashSet();
		for (Long kmer : path) {
			reads.addAll(kmers.get(kmer).getSupportingReads());
		}
		return reads;
	}
	private byte[] rescaleBaseQualities(List<Long> bases) {
		//Long largest = Collections.max(bases);
		//float scaleFactor = Math.min(1, MAX_QUAL_SCORE / (float)largest);
		byte[] result = new byte[bases.size()];
		for (int i = 0; i < result.length; i++) {
			//result[i] = (byte)(bases.get(i) * scaleFactor);
			result[i] = (byte)(bases.get(i) > MAX_QUAL_SCORE ? MAX_QUAL_SCORE : bases.get(i));
		}
		return result;
	}
	private LinkedList<Long> greedyTraverse(Long start) {
		LinkedList<Long> path = new LinkedList<Long>();
		Set<Long> visited = new HashSet<Long>();
		path.add(start);
		visited.add(start);
		for (Long node = greedyPrevState(start, visited); node != null; node = greedyPrevState(node, visited)) {
			path.addFirst(node);
			visited.add(node);
		}
		for (Long node = greedyNextState(start, visited); node != null; node = greedyNextState(node, visited)) {
			path.addLast(node);
			visited.add(node);
		}
		return path;
	}
	private Long greedyNextState(Long state, Set<Long> visited) {
		long best = -1;
		Long bestNode = null;
		for (Long next : KmerEncodingHelper.nextStates(state, k)) {
			DeBruijnNode node = kmers.get(next);
			if (node != null && node.getWeight() > best && !visited.contains(next)) {
				bestNode = next;
				best = node.getWeight();
			}
		}
		return bestNode; 
	}
	private Long greedyPrevState(Long state, Set<Long> visited) {
		long best = -1;
		Long bestNode = null;
		for (Long next : KmerEncodingHelper.prevStates(state, k)) {
			DeBruijnNode node = kmers.get(next);
			if (node != null && node.getWeight() > best && !visited.contains(next)) {
				bestNode = next;
				best = node.getWeight();
			}
		}
		return bestNode; 
	}
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("De Bruijn graph: k=%d, %d kmers, %d start kmers\n", k, kmers.size(), startkmers.size()));
		int max = 10;
		for (Long x : kmers.keySet()) {
			sb.append(String.format("%s(%d): %d weight from %d reads %s",
					KmerEncodingHelper.toString(k, x),
					x,
					kmers.get(x).getWeight(),
					kmers.get(x).getSupportingReads().size(),
					startkmers.containsKey(x) ? "(start)" : ""
					));
			sb.append(" from:{");
			for (Long y : KmerEncodingHelper.prevStates(x, k)) {
				DeBruijnNode node = kmers.get(y);
				if (node != null) {
					sb.append(KmerEncodingHelper.toString(k, y));
					sb.append(',');
				}
			}
			sb.append("} to:{");
			for (Long y : KmerEncodingHelper.nextStates(x, k)) {
				DeBruijnNode node = kmers.get(y);
				if (node != null) {
					sb.append(KmerEncodingHelper.toString(k, y));
					sb.append(',');
				}
			}
			sb.append("}\n");
			max--;
			if (max <= 0) break;
		}
		return sb.toString();
	}
}
