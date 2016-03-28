package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongList;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.ImmutableTriple;

import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;

/**
 * Corrects misassembly due to incorporation of read pair evidence at multiple positions
 * @author Daniel Cameron
 *
 */
class MisassemblyFixer {
	private final List<KmerPathSubnode> contig;
	/**
	 * path kmer offset -> { path node offset, pretransition kmers, posttransition kmers }
	 */
	private final NavigableMap<Integer, ImmutableTriple<Integer, LongList, LongList>> contigTransitionOffsetLookup;
	/**
	 * kmer -> { start, end, path kmer offset }
	 * unordered for now: ideally this should be an IntIntervalTreeMultimap
	 */
	private final Long2ObjectMap<List<ImmutableTriple<Integer, Integer, Integer>>> contigOffsetLookup;
	public MisassemblyFixer(Collection<KmerPathSubnode> contig) {
		this.contig = Lists.newArrayList(contig);
		this.contigTransitionOffsetLookup = createContigTransitionOffsetLookup(this.contig);
		this.contigOffsetLookup = createContigOffsetLookup(this.contig);
	}
	private static NavigableMap<Integer, ImmutableTriple<Integer, LongList, LongList>> createContigTransitionOffsetLookup(List<KmerPathSubnode> contig) {
		NavigableMap<Integer, ImmutableTriple<Integer, LongList, LongList>> lookup = new TreeMap<Integer, ImmutableTriple<Integer, LongList, LongList>>();
		int snoffset = 0;
		for (int i = 0; i < contig.size() - 1; i++) {
			KmerPathSubnode sn = contig.get(i);
			LongArrayList snendkmers = new LongArrayList();
			snendkmers.add(sn.lastKmer());
			for (int j = 0; j < sn.node().collapsedKmerOffsets().size(); j++) {
				int offset = sn.node().collapsedKmerOffsets().getInt(j);
				if (offset == sn.length() - 1) {
					snendkmers.add(sn.node().collapsedKmers().getLong(j));
				}
			}
			KmerPathSubnode snext = contig.get(i + 1);
			LongArrayList snextstartkmers = new LongArrayList();
			snextstartkmers.add(snext.firstKmer());
			for (int j = 0; j < snext.node().collapsedKmerOffsets().size(); j++) {
				int offset = snext.node().collapsedKmerOffsets().getInt(j);
				if (offset == 0) {
					snextstartkmers.add(snext.node().collapsedKmers().getLong(j));
				}
			}
			lookup.put(snoffset + sn.length() - 1, new ImmutableTriple<Integer, LongList, LongList>(i, snendkmers, snextstartkmers));
			snoffset += sn.length();
		}
		return lookup;
	}
	private static Long2ObjectMap<List<ImmutableTriple<Integer, Integer, Integer>>> createContigOffsetLookup(Collection<KmerPathSubnode> contig) {
		Long2ObjectMap<List<ImmutableTriple<Integer, Integer, Integer>>> contigOffsetLookup = new Long2ObjectOpenHashMap<List<ImmutableTriple<Integer, Integer, Integer>>>();
		int snoffset = 0;
		for (KmerPathSubnode sn : contig) {
			for (int i = 0; i < sn.length(); i++) {
				contigOffsetLookupAdd(contigOffsetLookup, snoffset + i, sn.node().kmer(i), sn.firstStart() + i, sn.firstEnd() + i);
			}
			for (int j = 0; j < sn.node().collapsedKmerOffsets().size(); j++) {
				int i = sn.node().collapsedKmerOffsets().getInt(j);
				contigOffsetLookupAdd(contigOffsetLookup, snoffset + i, sn.node().collapsedKmers().getLong(j), sn.firstStart() + i, sn.firstEnd() + i);
			}
			snoffset += sn.length();
		}
		return contigOffsetLookup;
	}
	private static void contigOffsetLookupAdd(Long2ObjectMap<List<ImmutableTriple<Integer, Integer, Integer>>> contigOffsetLookup, int offset, long kmer, int start, int end) {
		List<ImmutableTriple<Integer, Integer, Integer>> list = contigOffsetLookup.get(kmer);
		if (list == null) {
			list = new ArrayList<ImmutableTriple<Integer,Integer,Integer>>();
			contigOffsetLookup.put(kmer, list);
		}
		list.add(new ImmutableTriple<Integer, Integer, Integer>(start, end, offset));
	}
	/**
	 * Reassembles the given contig ensuring a valid traversal path
	 * @param contig
	 * @param evidence
	 */
	public List<KmerPathSubnode> correctMisassignedEvidence(Collection<KmerEvidence> evidence) {
		// we're only concerned about misassembly of read pairs
		// as a conservative approximation to full OLC reassembly
		// we naively left/right align everything and truncate if we have
		// zero support for a node transition.
		List<KmerPathSubnode> left = asLeftAligned(evidence);
		List<KmerPathSubnode> right = asRightAligned(evidence);
		boolean leftAnchored = left.get(0).prev().stream().anyMatch(sn -> sn.node().isReference());
		boolean rightAnchored = right.get(right.size() - 1).next().stream().anyMatch(sn -> sn.node().isReference());
		if (leftAnchored && !rightAnchored) return left;
		if (rightAnchored && !leftAnchored) return right;
		int leftWeight = left.stream().mapToInt(sn -> sn.weight()).sum();
		int rightWeight = right.stream().mapToInt(sn -> sn.weight()).sum();
		if (leftWeight >= rightWeight) return left;
		else return right;
	}
	private List<KmerPathSubnode> asLeftAligned(Collection<KmerEvidence> evidence) {
		int[] transitionSupport = transitionSupport(evidence, true);
		List<KmerPathSubnode> newContig = new ArrayList<KmerPathSubnode>(contig.size());
		newContig.add(contig.get(0));
		for (int i = 0; i < transitionSupport.length; i++) {
			if (transitionSupport[i] > 0) {
				newContig.add(contig.get(i + 1));
			} else {
				break;
			}
		}
		return newContig;
	}
	private List<KmerPathSubnode> asRightAligned(Collection<KmerEvidence> evidence) {
		int[] transitionSupport = transitionSupport(evidence, false);
		LinkedList<KmerPathSubnode> newContig = new LinkedList<KmerPathSubnode>();
		newContig.add(contig.get(contig.size() - 1));
		for (int i = contig.size() - 2; i >= 0; i--) {
			if (transitionSupport[i] > 0) {
				newContig.addFirst(contig.get(i));
			} else {
				break;
			}
		}
		return newContig;
	}
	private int[] transitionSupport(Collection<KmerEvidence> evidence, boolean leftAlign) {
		int[] transitionSupport = new int[contig.size() - 1];
		for (KmerEvidence e : evidence) {
			int[] offsets = matchingOffsets(e);
			int offset = leftAlign ? offsets[0] : offsets[1];
			addSupport(transitionSupport, e, offset);
		}
		return transitionSupport;
	}
	private void addSupport(int[] transitionSupport, KmerEvidence evidence, int offset) {
		for (Entry<Integer, ImmutableTriple<Integer, LongList, LongList>> entry : contigTransitionOffsetLookup.subMap(offset, offset + evidence.length() - 1).entrySet()) {
			long readKmerPreTransition = evidence.kmer(entry.getKey() - offset);
			long readKmerPostTransition = evidence.kmer(entry.getKey() - offset + 1);
			ImmutableTriple<Integer, LongList, LongList> transition = entry.getValue();
			if (transition.middle.contains(readKmerPreTransition) && transition.right.contains(readKmerPostTransition)) {
				transitionSupport[transition.left]++;
			}
		}
	}
	/**
	 * Returns the inferred contig offset of the starting read kmer for every kmer match  
	 * @param e read
	 * @return max and min best read starting kmer offset   
	 */
	private int[] matchingOffsets(KmerEvidence evidence) {
		SortedMultiset<Integer> counts = TreeMultiset.create();
		for (int i = 0; i < evidence.length(); i++) {
			KmerSupportNode n = evidence.node(i);
			if (n != null) {
				offsetLookup(counts, i, n.firstKmer(), n.firstStart(), n.firstEnd());
			}
		}
		int maxCount = counts.entrySet().stream()
				.mapToInt(e -> e.getCount())
				.max()
				.getAsInt();
		List<Integer> best = counts.entrySet().stream()
				.filter(e -> e.getCount() == maxCount)
				.map(e -> e.getElement())
				.collect(Collectors.toList());
		return new int[] { best.get(0), best.get(best.size() - 1) };
	}
	private void offsetLookup(Multiset<Integer> counts, int readOffset, long kmer, int positionStart, int positionEnd) {
		List<ImmutableTriple<Integer, Integer, Integer>> validPositionRanges = contigOffsetLookup.get(kmer);
		if (validPositionRanges != null) {
			for (ImmutableTriple<Integer, Integer, Integer> entry : validPositionRanges)  {
				if (IntervalUtil.overlapsClosed(entry.left, entry.middle, positionStart, positionEnd)) {
					counts.add(entry.right - readOffset);	
				}
			}
		}
	} 
}