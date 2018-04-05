package au.edu.wehi.idsv.debruijn.positional;

import java.util.BitSet;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;

import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.debruijn.ContigCategorySupportHelper;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;

public class PositionalContigCategorySupportHelper extends ContigCategorySupportHelper {
	public static String getCategorySupport(Collection<KmerPathSubnode> fullContig, Set<KmerEvidence> evidence, int k) {
		int kmersInPath = fullContig.stream()
				.mapToInt(x -> x.length())
				.sum();
		int categories = evidence.stream()
				.mapToInt(x -> ((SAMEvidenceSource)x.evidence().getEvidenceSource()).getSourceCategory())
				.max()
				.orElse(0) + 1;
		Long2ObjectOpenHashMap<RangeMap<Integer, Integer>> lookup = buildLookup(fullContig);
		List<BitSet> supportedKmers = Stream.generate(() -> new BitSet(kmersInPath))
				.limit(categories)
				.collect(Collectors.toList());
		for (KmerEvidence e : evidence) {
			int category = ((SAMEvidenceSource)e.evidence().getEvidenceSource()).getSourceCategory();
			BitSet support = supportedKmers.get(category);
			setSupportedBits(support, e, lookup);
		}
		List<BitSet> supportedBases = supportedKmers.stream()
			.map(x -> asSupportedBases(x, k))
			.collect(Collectors.toList());
		String cigars = asSupportCigars(kmersInPath + k - 1, supportedBases);
		return cigars;
	}
	private static void setSupportedBits(BitSet bs, KmerEvidence e, Long2ObjectOpenHashMap<RangeMap<Integer, Integer>> lookup) {
		for (int i = 0; i < e.length(); i++) {
			setSupportedBits(bs, e.kmer(i), e.startPosition() + i, e.endPosition() + i, lookup);
		}
	}
	private static void setSupportedBits(BitSet bs, long kmer, int start, int end, Long2ObjectOpenHashMap<RangeMap<Integer, Integer>> lookup) {
		RangeMap<Integer, Integer> rm = lookup.get(kmer);
		if (rm == null) return;
		RangeMap<Integer, Integer> srm = rm.subRangeMap(Range.closedOpen(start, end + 1));
		for (Integer offset : srm.asMapOfRanges().values()) {
			bs.set(offset);
		}
	}
	private static Long2ObjectOpenHashMap<RangeMap<Integer, Integer>> buildLookup(Collection<KmerPathSubnode> fullContig) {
		Long2ObjectOpenHashMap<RangeMap<Integer, Integer>> lookup = new Long2ObjectOpenHashMap<>(512);
		int offset = 0;
		for (KmerPathSubnode n : fullContig) {
			// primary kmers
			for (int i = 0; i < n.length(); i++) {
				addToLookup(offset + i, n.kmer(i), n.firstStart() + i, n.firstEnd() + i, lookup);
			}
			// error corrected kmers
			LongArrayList kmers = n.node().collapsedKmers();
			IntArrayList offsets = n.node().collapsedKmerOffsets();
			for (int i = 0; i < kmers.size(); i++) {
				addToLookup(offset + offsets.getInt(i), kmers.getLong(i), n.firstStart() + offsets.getInt(i), n.firstEnd() + offsets.getInt(i), lookup);
			}
			offset += n.length();
		}
		return lookup;
	}
	private static void addToLookup(int offset, long kmer, int start, int end, Long2ObjectOpenHashMap<RangeMap<Integer, Integer>> lookup) {
		RangeMap<Integer, Integer> rm = lookup.get(kmer);
		if (rm == null) {
			rm = TreeRangeMap.create();
			lookup.put(kmer, rm);
		}
		rm.put(Range.closedOpen(start, end + 1), offset);
	}
	// convert per-kmer to per base: kmer only support intra-kmer base transitions
	private static BitSet asSupportedBases(BitSet supportedKmers, int k) {
		BitSet bases = new BitSet(supportedKmers.size() + k - 1);
		for (int i = 0; i < supportedKmers.size(); i++) {
			if (supportedKmers.get(i)) {
				// example, k = 3
				// 0123456
				// ^
				// at i=0, breakpoints before 1 & 2 are supported
				bases.set(i + 1, i + k);
			}
		}
		return bases;
	}
}
