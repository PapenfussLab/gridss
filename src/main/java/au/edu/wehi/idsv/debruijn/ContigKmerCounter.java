package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
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
import java.util.stream.Stream;

/**
 * Counts matching kmers
 */
public class ContigKmerCounter {
    private static final Log log = Log.getInstance(ContigKmerCounter.class);
    private final List<String> contigs = new ArrayList<>();
    private final KmerLookup kmerLookup;
    private final LongList counts = new LongArrayList();
    private final int k;
    public ContigKmerCounter(List<String> contigs, List<byte[]> sequences, int k, int stride) {
        this(deduplicatedReference(contigs, sequences), k, stride);
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
    public ContigKmerCounter(List<ReferenceSequence> ref, int k, int stride) {
        this(
            ref.stream().map(rs -> rs.getName()).collect(Collectors.toList()),
            ref.stream().map(rs -> rs.getBases()).collect(Collectors.toList()),
            k,
            stride);
    }
    public ContigKmerCounter(TwoBitBufferedReferenceSequenceFile reference, int k, int stride) {
        this.k = k;
        this.kmerLookup = k <= 16 ? new KmerLookup32() : new KmerLookup64();
        for (SAMSequenceRecord ssr : reference.getSequenceDictionary().getSequences()) {
            String contig = ssr.getContig();
            TwoBitBufferedReferenceSequenceFile.PackedReferenceSequence prs = reference.getPackedSequence(contig);
            log.debug("Adding:\t" + contig);
            sequentialAddToLookup(contig, prs, stride);
        }
    }
    public ContigKmerCounter(Stream<ReferenceSequence> ref, int k, int stride) {
        this.k = k;
        this.kmerLookup = k <= 16 ? new KmerLookup32() : new KmerLookup64();
        ref.forEachOrdered(rs -> {
            sequentialAddToLookup(rs.getName(), new TwoBitBufferedReferenceSequenceFile.PackedReferenceSequence(rs), stride);
        });
    }

    private void sequentialAddToLookup(String contig, TwoBitBufferedReferenceSequenceFile.PackedReferenceSequence prs, int stride) {
        if (contigs.contains(contig)) return; // prevent double-counting
        contigs.add(contig);
        counts.add(0);
        int contigId = contigs.size() - 1;
        for (int i = 0; i < prs.length() - k + 1; i += stride) {
            int startRefPos = i + 1;
            int endRefPos = i + k;
            if (!prs.anyAmbiguous(startRefPos, endRefPos)) {
                long kmer = prs.getKmer(i, k);
                sequentialAddToLookup(contigId, kmer);
                sequentialAddToLookup(contigId, KmerEncodingHelper.reverseComplement(k, kmer));
            }
        }
    }

    private void sequentialAddToLookup(int contigId, long kmer) {
        IntList value = kmerLookup.getKmer(kmer);
        if (value == null) {
            value = new IntArrayList(4);
            value.add(contigId);
            kmerLookup.putKmer(kmer, value);
        } else if (value.getInt(value.size() - 1) != contigId) { // repeated kmers in the reference shouldn't count multiple times
            value.add(contigId);
        }
    }

    public int count(byte[] seq) {
        int hits = 0;
        PackedSequence fps = new PackedSequence(seq, false, false);
        //PackedSequence bps = new PackedSequence(seq, true, true); // don't need to RC the reads since we added the RC of the reference
        for (int i = 0; i < seq.length - k + 1; i++) {
            hits += count(fps.getKmer(i, k));
            //hits += count(bps.getKmer(i, k));
        }
        return hits;
    }

    private int count(long kmer) {
        IntList hits = kmerLookup.getKmer(kmer);
        if (hits == null) return 0;
        hits.forEach((int refIndex) -> counts.set(refIndex, counts.getLong(refIndex) + 1));
        return hits.size();
    }

    public List<String> getContigs() {
        return contigs;
    }
    public LongList getKmerCounts() {
        return counts;
    }
    private abstract static class KmerLookup {
        protected abstract void putKmer(long kmer, IntList list);

        protected abstract IntList getKmer(long kmer);
    }
    private static class KmerLookup32 extends KmerLookup {
        private final Int2ObjectOpenHashMap<IntList> kmerLookup = new Int2ObjectOpenHashMap<>();
        @Override
        protected void putKmer(long kmer, IntList list) {
            kmerLookup.put((int)kmer, list);
        }

        @Override
        protected IntList getKmer(long kmer) {
            return kmerLookup.get((int)kmer);
        }
    }
    private static class KmerLookup64 extends KmerLookup {
        private final Long2ObjectOpenHashMap<IntList> kmerLookup = new Long2ObjectOpenHashMap<>();
        @Override
        protected void putKmer(long kmer, IntList list) {
            kmerLookup.put(kmer, list);
        }

        @Override
        protected IntList getKmer(long kmer) { return kmerLookup.get(kmer); }
    }
}
