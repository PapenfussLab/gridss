package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongList;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Counts matching kmers
 */
public class ContigKmerCounter {
    private final List<String> contigs = new ArrayList<>();
    private final Long2ObjectOpenHashMap<IntList> kmerLookup = new Long2ObjectOpenHashMap<>();
    private final LongList counts = new LongArrayList();
    private final int k;
    public ContigKmerCounter(List<String> contigs, List<byte[]> sequences, int k) {
        this(deduplicatedReference(contigs, sequences), k);
    }
    private static TwoBitBufferedReferenceSequenceFile deduplicatedReference(List<String> contigs, List<byte[]> sequences) {
        contigs = Lists.newArrayList(contigs);
        sequences = Lists.newArrayList(sequences);
        Set<String> found = new HashSet<>();
        for (int i = 0; i < contigs.size(); i++) {
            if (found.contains(contigs.get(i))) {
                contigs.remove(i);
                sequences.remove(i);
                i--;
            } else {
                found.add(contigs.get(i));
            }
        }
        return new TwoBitBufferedReferenceSequenceFile(new InMemoryReferenceSequenceFile(contigs, sequences));
    }
    public ContigKmerCounter(List<ReferenceSequence> ref, int k) {
        this(
            ref.stream().map(rs -> rs.getName()).collect(Collectors.toList()),
            ref.stream().map(rs -> rs.getBases()).collect(Collectors.toList()),
            k);
    }
    public ContigKmerCounter(TwoBitBufferedReferenceSequenceFile reference, int k) {
        this.k = k;
        for (SAMSequenceRecord ssr : reference.getSequenceDictionary().getSequences()) {
            String contig = ssr.getContig();
            TwoBitBufferedReferenceSequenceFile.PackedReferenceSequence prs = reference.getPackedSequence(contig);
            addToLookup(contig, prs);
        }
    }

    private void addToLookup(String contig, TwoBitBufferedReferenceSequenceFile.PackedReferenceSequence prs) {
        if (contigs.contains(contig)) return; // prevent double-counting
        contigs.add(contig);
        counts.add(0);
        int contigId = contigs.size() - 1;
        for (int i = 0; i < prs.length() - k + 1; i++) {
            int startRefPos = i + 1;
            int endRefPos = i + k;
            if (!prs.anyAmbiguous(startRefPos, endRefPos)) {
                long kmer = prs.getKmer(i, k);
                addToLookup(contigId, kmer);
                addToLookup(contigId, KmerEncodingHelper.reverseComplement(k, kmer));
            }
        }
    }

    private void addToLookup(int contigId, long kmer) {
        IntList value = kmerLookup.get(kmer);
        if (value == null) {
            value = new IntArrayList(2);
            kmerLookup.put(kmer, value);
        }
        if (!value.contains(contigId)) {
            // repeated kmers in the reference shouldn't count multiple times
            // we're only interested in inclusion
            value.add(contigId);
        }
    }

    public int count(byte[] seq) {
        int hits = 0;
        PackedSequence fps = new PackedSequence(seq, false, false);
        //PackedSequence bps = new PackedSequence(seq, true, true); // don't need to RC the reads since we add the RC of the reference
        for (int i = 0; i < seq.length - k + 1; i++) {
            hits += count(fps.getKmer(i, k));
            //hits += count(bps.getKmer(i, k));
        }
        return hits;
    }

    private int count(long kmer) {
        IntList hits = kmerLookup.get(kmer);
        if (hits == null) return 0;
        for (int i = 0; i < hits.size(); i++) {
            counts.set(hits.getInt(i), counts.getLong(i) + 1);
        }
        return hits.size();
    }

    public List<String> getContigs() {
        return contigs;
    }
    public LongList getKmerCounts() {
        return counts;
    }
}
