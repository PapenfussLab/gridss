package au.edu.wehi.idsv.picard;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Log;

import java.io.IOException;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.PackedSequence;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

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
	private volatile ImmutableMap<String, PackedReferenceSequence> cache = ImmutableMap.of();
	private PackedReferenceSequence[] referenceIndexLookup; 
	
	public TwoBitBufferedReferenceSequenceFile(ReferenceSequenceFile underlying) {
		this.underlying = underlying;
		this.referenceIndexLookup = new PackedReferenceSequence[underlying.getSequenceDictionary().getSequences().size()];
	}
	public byte getBase(int referenceIndex, int position) {
		PackedReferenceSequence seq = referenceIndexLookup[referenceIndex];
		if (seq == null) {
			synchronized (referenceIndexLookup) {
				seq = addToCache(underlying.getSequenceDictionary().getSequence(referenceIndex).getSequenceName());
			}
		}
		return seq.get(position - 1); 
	}
	private class PackedReferenceSequence extends PackedSequence {
		private final String name;
	    private final int contigIndex;
	    private final long length;
	    /**
	     * 1-based lookup of ambiguous bases
	     */
	    private final RangeSet<Integer> ambiguous = TreeRangeSet.create();
		public PackedReferenceSequence(ReferenceSequence seq) {
			super(seq.getBases(), false, false);
			this.name = seq.getName();
			this.contigIndex = seq.getContigIndex();
			this.length = seq.length();
			for (int i = 0; i < length; i++) {
				if (KmerEncodingHelper.isAmbiguous(seq.getBases()[i])) {
					int end = i + 1;
					while (end < length && KmerEncodingHelper.isAmbiguous(seq.getBases()[end])) {
						end++;
					}
					ambiguous.add(Range.closedOpen(i + 1, end + 1));
					i = end - 1;
				}
			}
		}
		public ReferenceSequence getSequence() {
			return getSubsequenceAt(1, length);
		}
		public ReferenceSequence getSubsequenceAt(long start, long stop) {
			int length = (int)(stop - start + 1);
			ReferenceSequence seq = new ReferenceSequence(name, contigIndex, getBytes((int)(start - 1), length));
			for (Range<Integer> r : ambiguous.subRangeSet(Range.closedOpen((int)start, (int)start + length)).asRanges()) {
				for (int i = r.lowerEndpoint(); i < r.upperEndpoint(); i++) {
					int offset = (int)(i - start);
					seq.getBases()[offset] = 'N';
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
	private synchronized PackedReferenceSequence addToCache(String contig) {
		PackedReferenceSequence seq = cache.get(contig);
		if (seq != null) {
			// already populated by another thread while we were waiting to enter
			// this synchronized block
			return seq;
		}
		log.debug("Caching reference genome contig ", contig);
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
