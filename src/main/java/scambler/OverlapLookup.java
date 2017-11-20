package scambler;

import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.debruijn.PackedSequence;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

public class OverlapLookup {
	private final int MAX_MISMATCHES = 0;
	private final int minOverlap;
	private final int kmerSize;
	private final Long2ObjectOpenHashMap<List<Read>> lookup = new Long2ObjectOpenHashMap<>();
	public OverlapLookup(int minOverlap) {
		this.minOverlap = minOverlap;
		this.kmerSize = Math.min(32, minOverlap);
	}
	public void add(Read r) {
		long key = r.getSeq().getKmer(0, kmerSize);
		List<Read> entry = lookup.get(key);
		if (entry == null) {
			entry = new ArrayList<>(4);
			lookup.put(key, entry);
		}
		entry.add(r);
	}
	public List<Overlap> successors(Read r) {
		PackedSequence seq = r.getSeq();
		List<Overlap> overlaps = new ArrayList<>();
		for (int i = 0; i < seq.length() - kmerSize; i++) {
			long kmer = seq.getKmer(i, kmerSize);
			for (Read hit : lookup.get(kmer)) {
				if (hit != r) {
					Overlap o = new Overlap(r, hit, i);
					if (o.matchingBases + MAX_MISMATCHES >= o.overlap && o.overlap >= minOverlap) {
						overlaps.add(o);
					}
				}
			}
		}
		return overlaps;
	}
}
