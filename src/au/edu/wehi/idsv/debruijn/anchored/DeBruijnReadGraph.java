package au.edu.wehi.idsv.debruijn.anchored;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.ReadKmerIterable;
import au.edu.wehi.idsv.sam.AnomolousReadAssembly;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;

public class DeBruijnReadGraph extends DeBruijnGraphBase<DeBruijnNodeBase> {
	private final Multimap<Long, Integer> startkmers = HashMultimap.<Long, Integer>create();
	public DeBruijnReadGraph(int k, BreakendDirection direction) {
		super(k, direction);
	}
	@Override
	protected DeBruijnNodeBase createEmptyNode() {
		return new DeBruijnNodeBase();
	}
	@Override
	protected void addKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super.addKmer(evidence, readKmerOffset, kmer);
		if (evidence.getReferenceKmerCount() > 0) {
			// we have a fully reference anchored kmer
			if (evidence.isReferenceAnchor(readKmerOffset)) {
				startkmers.put(kmer.kmer, 0);
			}
		} else {
			// no kmer is fully anchored, use the first kmer as that has the most anchored bases
			if (evidence.isFirstKmer(readKmerOffset) && evidence.basesSupportingReference(readKmerOffset) > 0) {
				int offset = k - evidence.basesSupportingReference(readKmerOffset);
				startkmers.put(kmer.kmer, offset);
			}
		}
	}
	@Override
	protected void removeKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		super.removeKmer(evidence, readKmerOffset, kmer);
		if (evidence.getReferenceKmerCount() > 0) {
			// we have a fully reference anchored kmer
			if (evidence.isReferenceAnchor(readKmerOffset)) {
				startkmers.remove(kmer.kmer, 0);
			}
		} else {
			// no kmer is fully anchored, use the first kmer as that has the most anchored bases
			if (evidence.isFirstKmer(readKmerOffset) && evidence.basesSupportingReference(readKmerOffset) > 0) {
				int offset = k - evidence.basesSupportingReference(readKmerOffset);
				startkmers.remove(kmer.kmer, offset);
			}
		}
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
			int offsetSize = startkmers.get(breakpointAnchor).iterator().next();
			softclipSize += offsetSize;
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
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(super.toString());
		sb.append(String.format("%d start kmers", startkmers.size()));
		return sb.toString();
	}
}
