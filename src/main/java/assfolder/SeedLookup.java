package assfolder;

import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

import java.util.*;

public class SeedLookup {
    private final int k;
    private final Long2ObjectOpenHashMap<Collection<ReadOffset>> lookup = new Long2ObjectOpenHashMap<>();

    public SeedLookup(int k) {
        if (k > 32) throw new IllegalArgumentException("kmer seed size must fit in 64-bit integer");
        this.k = k;
    }

    public Collection<ReadOffset> findOverlaps(Read read) {
        Set<ReadOffset> matches = new HashSet<>();
        for (int i = 0 ; i < read.length() - k + 1; i++) {
            long kmer = read.getKmer(i, k);
            Collection<ReadOffset> coll = lookup.get(kmer);
            if (coll != null) {
                for (ReadOffset ro : coll) {
                    if (ro.read != read) {
                        ReadOffset hit = new ReadOffset(ro.read, i - ro.offset);
                        matches.add(hit);
                    }
                }
            }
        }
        return matches;
    }

    public void add(Read read) {
        for (int i = 0 ; i < read.length() - k + 1; i++) {
            ReadOffset ro = new ReadOffset(read, i);
            long kmer = read.getKmer(i, k);
            Collection<ReadOffset> coll = lookup.get(kmer);
            if (coll == null) {
                coll = new ArrayList<>();
                lookup.put(kmer, coll);
            }
            coll.add(ro);
        }
    }

    public void remove(Read read) {
        for (int i = 0 ; i < read.length() - k + 1; i++) {
            long kmer = read.getKmer(i, k);
            Collection<ReadOffset> coll = lookup.get(kmer);
            if (coll == null) {
                throw new NoSuchElementException("Unable to find read kmer in lookup");
            }
            coll.removeIf(ro -> ro.read == read);
            if (coll.isEmpty()) {
                lookup.remove(kmer);
            }
        }
    }

    /**
     * Finds the read in the graph (if it already exists).
     * @param read
     * @return
     */
    public Read getRead(Read read) {
        long kmer = read.getKmer(0, k);
        Collection<ReadOffset> coll = lookup.get(kmer);
        if (coll == null) return null;
        for (ReadOffset ro : coll) {
            if (ro.read.uid().equals(read.uid())) {
                return ro.read;
            }
        }
        return null;
    }
}
