package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.debruijn.ContigCategorySupportHelper;
import com.google.common.collect.*;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import org.apache.commons.lang3.tuple.Pair;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class PositionalContigCategorySupportHelper extends ContigCategorySupportHelper {
	/**
	 * Calculates the level of support for a breakpoint at the given contig position.
	 * Split read are considered to support all positions in the contig that the read overlaps.
	 * Read pairs are considered to support all position between the start of the contig and read contig position.
	 * @param fullContig full contig path
	 * @param evidence supporting evidence
	 * @param k kmer size
	 * @param assemblyDirection direction of assembly
	 * @return Returns a pair of per category support lists.
	 * Each entry in the first list contain the number of supporting reads for each breakpoint position.
	 * Each entry in the second list contain the corresponding net quality scores for each breakpoint position.
	 */
	public static Pair<List<int[]>, List<float[]>> getPerCategoryPerBaseSupport(Collection<KmerPathSubnode> fullContig, int categories, Stream<KmerEvidence> evidence, int k, BreakendDirection assemblyDirection, AssemblyConfiguration ac) {
		int kmersInPath = fullContig.stream()
				.mapToInt(x -> x.length())
				.sum();
		Long2ObjectOpenHashMap<RangeMap<Integer, Integer>> lookup = buildLookup(fullContig);

		// Counts and qual of evidence supporting a breakpoint immediately before the base at the given position
		List<int[]> supportingEvidenceCount = Stream.generate(() -> new int[kmersInPath + k])
				.limit(categories)
				.collect(Collectors.toList());
		List<float[]> supportingEvidenceQual = Stream.generate(() -> new float[kmersInPath + k])
				.limit(categories)
				.collect(Collectors.toList());
		evidence.forEach(e -> {
			int cat = ((SAMEvidenceSource)e.evidence().getEvidenceSource()).getSourceCategory();
			if (e.evidence() instanceof SingleReadEvidence) {
				trackSingleReadEvidence(k, lookup, supportingEvidenceCount.get(cat), supportingEvidenceQual.get(cat), e);
			} else {
				assert(e.evidence() instanceof NonReferenceReadPair);
				trackReadPairEvidence(k, lookup, supportingEvidenceCount.get(cat), supportingEvidenceQual.get(cat), e, assemblyDirection, ac);
			}
		});
		return Pair.of(supportingEvidenceCount, supportingEvidenceQual);
	}

	private static void trackReadPairEvidence(
			int k,
			Long2ObjectOpenHashMap<RangeMap<Integer,Integer>> lookup,
			int[] supportingEvidenceCount,
			float[] supportingEvidenceQual,
			KmerEvidence e,
			BreakendDirection assemblyDirection,
			AssemblyConfiguration ac) {
		Range<Integer> bounds = contigBaseOffsetBounds(lookup, e);
		Range<Integer> anchorBounds = null;
		if (ac.includePairAnchors) {
			NonReferenceReadPair nrrp = (NonReferenceReadPair)e.evidence();
			KmerEvidence e2 = KmerEvidence.createAnchor(k, nrrp, ac.pairAnchorMismatchIgnoreEndBases, nrrp.getEvidenceSource().getContext().getReference());
			anchorBounds = contigBaseOffsetBounds(lookup, e2);
		}
		if (anchorBounds == null) {
			if (assemblyDirection == BreakendDirection.Forward) {
				bounds = Range.closed(0, bounds.upperEndpoint());
			} else {
				bounds = Range.closed(bounds.lowerEndpoint(), supportingEvidenceCount.length);
			}
		} else {
			bounds = Range.closed(Math.min(bounds.lowerEndpoint(), anchorBounds.lowerEndpoint()), Math.min(bounds.upperEndpoint(), anchorBounds.upperEndpoint()));
		}
		for (int i = bounds.lowerEndpoint(); i <= bounds.upperEndpoint() + k - 1; i++) {
			supportingEvidenceCount[i]++;
			supportingEvidenceQual[i] += e.evidenceQuality();
		}
	}
	private static Range<Integer> contigBaseOffsetBounds(Long2ObjectOpenHashMap<RangeMap<Integer,Integer>> lookup, KmerEvidence e) {
		int min = Integer.MAX_VALUE;
		int max = Integer.MIN_VALUE;
		for (int i = 0; i < e.length(); i++) {
			for (Integer contigKmerOffset : getContigBaseOffsetFor(lookup, e, i)) {
				min = Math.min(contigKmerOffset, min);
				max = Math.max(contigKmerOffset, max);
			}
		}
		if (min == Integer.MAX_VALUE) {
			return null;
		}
		return Range.closed(min, max);
	}
	private static void trackSingleReadEvidence(
			int k,
			Long2ObjectOpenHashMap<RangeMap<Integer,Integer>> lookup,
			int[] supportingEvidenceCount,
			float[] supportingEvidenceQual,
			KmerEvidence e) {
		RangeSet<Integer> supportedBaseOffsets = TreeRangeSet.create();
		for (int i = 0; i < e.length(); i++) {
			for (Integer contigKmerOffset : getContigBaseOffsetFor(lookup, e, i)) {
				// read kmer support is support for just the base transitions within the kmer
				supportedBaseOffsets.add(Range.closedOpen(contigKmerOffset + 1, contigKmerOffset + k));
			}
		}
		for (Range<Integer> r : supportedBaseOffsets.asRanges()) {
			for (int i = r.lowerEndpoint(); i < r.upperEndpoint(); i++) {
				supportingEvidenceCount[i]++;
				supportingEvidenceQual[i] += e.evidenceQuality();
			}
		}
	}

	/**
	 * Determins where offset in the given evidence is included in the assembly.
	 * Read pairs can overlap multiple times if the sequence kmer is repeated.
	 * Since we don't actually know which position we placed the read in, we'll return them all.
	 */
	private static Collection<Integer> getContigBaseOffsetFor(Long2ObjectOpenHashMap<RangeMap<Integer,Integer>> lookup, KmerEvidence e, int evidenceOffset) {
		KmerSupportNode node = e.node(evidenceOffset);
		if (node != null) {
			long kmer = node.firstKmer();
			RangeMap<Integer, Integer> lookupEntry = lookup.get(kmer);
			if (lookupEntry != null) {
				Map<Range<Integer>, Integer> matching = lookupEntry.subRangeMap(Range.closedOpen(node.firstStart(), node.firstEnd() + 1)).asMapOfRanges();
				if (!matching.isEmpty()) {
					return matching.values();
				}
			}
		}
		return Collections.emptyList();
	}

	private static int[] kmerToBaseSupport(int k, int[] a) {
		int[] result = new int[a.length + k - 1];
		System.arraycopy(a, 0, result, 0, a.length);
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < k; j++) {
				result[i + j] = Math.max(result[i], result[i + j]);
			}
		}
		return result;
	}
	private static float[] kmerToBaseSupport(int k, float[] a) {
		float[] result = new float[a.length + k - 1];
		System.arraycopy(a, 0, result, 0, a.length);
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < k; j++) {
				result[i + j] = Math.max(result[i], result[i + j]);
			}
		}
		return result;
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
}
