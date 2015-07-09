package au.edu.wehi.idsv.picard;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Log;

import java.io.IOException;
import java.util.BitSet;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.PackedSequence;

import com.google.common.collect.ImmutableMap;

/**
 * 2bit encodes and buffers the entire reference to enable efficient random lookup of small subsequences
 * @author cameron.d
 *
 */
public class TwoBitBufferedReferenceSequenceFile implements ReferenceSequenceFile, ReferenceLookup {
	private static final Log log = Log.getInstance(TwoBitBufferedReferenceSequenceFile.class);
	private final ReferenceSequenceFile underlying;
	/**
	 * Cached contigs
	 */
	private ImmutableMap<String, PackedReferenceSequence> cache = ImmutableMap.of();
	private PackedReferenceSequence[] referenceIndexLookup;
	public TwoBitBufferedReferenceSequenceFile(ReferenceSequenceFile underlying) {
		this.underlying = underlying;
		this.referenceIndexLookup = new PackedReferenceSequence[underlying.getSequenceDictionary().getSequences().size()];
		log.debug("Loading reference genome into memory");
		for (SAMSequenceRecord contig : underlying.getSequenceDictionary().getSequences()) {
			getSequence(contig.getSequenceName());
		}
		log.debug("Reference genome loaded");
	}
	public byte getBase(int referenceIndex, int position) {
		PackedReferenceSequence seq = referenceIndexLookup[referenceIndex];
		if (seq.ambiguous.get(position - 1)) {
			return 'N';
		}
		return seq.get(position - 1);
	}
	private class PackedReferenceSequence extends PackedSequence {
		private final String name;
	    private final int contigIndex;
	    private final long length;
	    private final BitSet ambiguous; 
		public PackedReferenceSequence(ReferenceSequence seq) {
			super(seq.getBases(), false, false);
			this.name = seq.getName();
			this.contigIndex = seq.getContigIndex();
			this.length = seq.length();
			this.ambiguous = new BitSet(seq.length());
			byte[] seqBases = seq.getBases();
			for (int i = 0; i < length; i++) {
				if (KmerEncodingHelper.isAmbiguous(seqBases[i])) {
					ambiguous.set(i);
				}
			}
		}
		public ReferenceSequence getSequence() {
			return getSubsequenceAt(1, length);
		}
		public ReferenceSequence getSubsequenceAt(long start, long stop) {
			int length = (int)(stop - start + 1);
			ReferenceSequence seq = new ReferenceSequence(name, contigIndex, getBytes((int)(start - 1), length));
			byte[] seqBases = seq.getBases();
			for (int i = 0; i < length; i++) {
				if (ambiguous.get((int)start - 1 + i)) {
					seqBases[i] = 'N';
				}
			}
			return seq;
		}
	}
	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		return underlying.getSequenceDictionary();
	}
	@Override
	public ReferenceSequence nextSequence() {
		return underlying.nextSequence();
	}
	@Override
	public void reset() {
		underlying.reset();
	}
	@Override
	public boolean isIndexed() {
		return underlying.isIndexed();
	}
	/**
	 * Updates the cache to include the new contig
	 * @param contig
	 */
	private PackedReferenceSequence addToCache(String contig) {
		PackedReferenceSequence seq = cache.get(contig);
		if (seq != null) {
			// already populated by another thread while we were waiting to enter
			// this synchronized block
			return seq;
		}
		ReferenceSequence fullContigSequence = underlying.getSequence(contig);
		seq = new PackedReferenceSequence(fullContigSequence);
		cache = ImmutableMap.<String, PackedReferenceSequence>builder()
				.putAll(cache)
				.put(contig, seq)
				.build();
		referenceIndexLookup[underlying.getSequenceDictionary().getSequence(contig).getSequenceIndex()] = seq;
		return seq;
	}
	@Override
	public ReferenceSequence getSequence(String contig) {
		PackedReferenceSequence seq = cache.get(contig);
		if (seq == null) {
			seq = addToCache(contig);
		}
		return seq.getSequence();
	}
	@Override
	public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
        PackedReferenceSequence seq = cache.get(contig);
		if (seq == null) {
			seq = addToCache(contig);
		}
		return seq.getSubsequenceAt(start, stop);
	}
	@Override
	public void close() throws IOException {
		underlying.close();	
	}
}
