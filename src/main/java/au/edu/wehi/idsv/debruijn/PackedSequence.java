package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.util.IntervalUtil;
import com.google.common.math.IntMath;

import java.io.Serializable;
import java.math.RoundingMode;

/**
 * Compresses the given sequence by representing in 2-bit format
 * @author Daniel Cameron
 *
 */
public class PackedSequence implements Serializable {
	private static final int BITS_PER_BASE = 2;
	private static final int BASES_PER_WORD = Long.SIZE / BITS_PER_BASE;
	private static final int ARRAY_SHIFT = Long.SIZE - 1 - Long.numberOfLeadingZeros(BASES_PER_WORD);
	private static final int ARRAY_OFFSET_MASK = (1 << ARRAY_SHIFT) - 1;
	private static final long BASE_MASK = (1 << BITS_PER_BASE) - 1;
	/**
	 * First base is packed in MSB of first word
	 * Second base is packed in second MSB of first word
	 * and so on
	 */
	private final long[] packed;
	private final int baseCount;
	public PackedSequence(long[] packed, int baseCount) {
		this.packed = packed;
		this.baseCount = baseCount;
	}
	public PackedSequence(byte[] bases, boolean reverse, boolean complement) {
		packed = new long[IntMath.divide(bases.length, BASES_PER_WORD, RoundingMode.CEILING)];
		baseCount = bases.length;
		int revMultiplier = 1;
		int revOffset = 0;
		if (reverse) {
			revMultiplier = -1;
			revOffset = bases.length - 1;
		}
		if (complement) { 
			for (int i = 0; i < bases.length; i++) {
				byte base = bases[i];
				long encoded = KmerEncodingHelper.picardBaseToEncoded(base);
				encoded = KmerEncodingHelper.complement(1, encoded);
				setBaseEncoded(revOffset + i * revMultiplier, encoded);
			}
		} else {
			for (int i = 0; i < bases.length; i++) {
				byte base = bases[i];
				long encoded = KmerEncodingHelper.picardBaseToEncoded(base);
				setBaseEncoded(revOffset + i * revMultiplier, encoded);
			}
		}
	}
	/**
	 * Subsequence of the given packed sequence
	 * @param offset 0-based position to start
	 * @param length length
	 */
	public PackedSequence(PackedSequence sequence, int offset, int length) {
		if (offset + length > sequence.length()) {
			throw new IllegalArgumentException("subsequence out of bounds");
		}
		this.packed = new long[IntMath.divide(length, BASES_PER_WORD, RoundingMode.CEILING)];
		this.baseCount = length;
		if (length == 0) { 
			return;
		}
		for (int i = 0; i < packed.length; i++) {
			int wordOffset = offset + i * BASES_PER_WORD;
			int remainingLength = length - i * BASES_PER_WORD;
			int wordLength = Math.min(BASES_PER_WORD, remainingLength);
			// bit shift is required to pack the sequence bases for the final word in the MSBs.
			this.packed[i] = sequence.getKmer(wordOffset, wordLength) << 2 * (BASES_PER_WORD - wordLength);
		}
	}
	/**
	 * Concatenates the given sequences
	 * @param seq1
	 * @param seq2
	 */
	public PackedSequence(PackedSequence seq1, PackedSequence seq2) {
		this.baseCount = seq1.baseCount + seq2.baseCount;
		this.packed = new long[IntMath.divide(this.baseCount, BASES_PER_WORD, RoundingMode.CEILING)];
		System.arraycopy(seq1.packed, 0, this.packed, 0, seq1.packed.length);
		int wordOffset = seq1.length() % BASES_PER_WORD;
		int wordIndex = seq1.length() / BASES_PER_WORD;
		if (wordOffset == 0) {
			// we can just copy seq2 since it's word aligned
			System.arraycopy(seq2.packed, 0, this.packed, wordIndex, seq2.packed.length);
		}
		for (int i = 0; i < seq2.packed.length; i++) {
			long word = seq2.packed[i];
			int thisWordBases = BASES_PER_WORD - wordOffset;
			int nextWordBases = wordOffset;
			long msb = word >>> (2 * nextWordBases);
			long lsb = word & ((1L << (2 * nextWordBases)) - 1);
			this.packed[wordIndex + i] |= msb;
			if (wordIndex + i< this.packed.length - 1) {
				this.packed[wordIndex + i + 1] |= lsb << (2 * thisWordBases);
			}
		}
	}
	private void setBaseEncoded(final int offset, final long base) {
		if (offset < 0 || offset >= baseCount) throw new IllegalArgumentException("offset must fall within sequence");
		int wordIndex = offset >> ARRAY_SHIFT;
		int wordOffset = BASES_PER_WORD - 1 - (offset & ARRAY_OFFSET_MASK);
		//assert(wordIndex < packed.length);
		long word = packed[wordIndex];
		word &= ~(BASE_MASK << (BITS_PER_BASE * wordOffset));
		word |= base << (BITS_PER_BASE * wordOffset);
		packed[wordIndex] = word;
	}
	private long getBaseEncoded(final int offset) {
		if (offset < 0 || offset >= baseCount) throw new IllegalArgumentException("offset must fall within sequence");
		int wordIndex = offset >> ARRAY_SHIFT;
		int wordOffset = BASES_PER_WORD - 1 - (offset & ARRAY_OFFSET_MASK);
		//assert(wordIndex < packed.length);
		long base = packed[wordIndex] >>> (BITS_PER_BASE * wordOffset);
		return base;
	}
	private long getWordBases(final int wordIndex, final int highBaseIgnoreCount, final int lowBaseIgnoreCount) {
		long word = packed[wordIndex];
		word <<= BITS_PER_BASE * highBaseIgnoreCount; // force high bases off the top
		word >>>= BITS_PER_BASE * (highBaseIgnoreCount + lowBaseIgnoreCount); // and low off the bottom
		return word;
	}
	public byte get(final int offset) {
		byte b = KmerEncodingHelper.encodedToPicardBase(getBaseEncoded(offset));
		return b;
	}
	public byte[] getBytes(int offset, int length) {
		assert(offset + length <= packed.length * BASES_PER_WORD);
		byte[] seq = new byte[length];
		for (int i = 0; i < length; i++) {
			seq[i] = get(offset + i);
		}
		return seq;
	}
	public long getKmer(final int offset, final int length) {
		assert(offset + length <= packed.length * BASES_PER_WORD);
		int wordIndex = offset >> ARRAY_SHIFT;
		int basesToSkipInWord = offset & ARRAY_OFFSET_MASK;
		int basesRemaining = BASES_PER_WORD - basesToSkipInWord;
		if (length <= BASES_PER_WORD - basesToSkipInWord) {
			return getWordBases(wordIndex, basesToSkipInWord, BASES_PER_WORD - basesToSkipInWord - length);
		} else {
			int lengthInNextWord = length - basesRemaining;
			long kmer = getWordBases(wordIndex, basesToSkipInWord, 0);
			kmer <<= lengthInNextWord * BITS_PER_BASE;
			kmer |= getWordBases(wordIndex + 1, 0, BASES_PER_WORD - lengthInNextWord);
			return kmer;
		}
	}
	@Override
	public String toString() {
		return new String(getBytes(0, baseCount));
	}
	private static final long MATCH_MASK = 0x5555555555555555L;
	public static int overlapMatches(PackedSequence s1, PackedSequence s2, int start2RelativeToStart1) {
		int matches = 0;
		int overlapLength = overlapLength(s1, s2, start2RelativeToStart1);
		if (overlapLength <= 0) {
			return 0;
		}
		// find common starting point
		int s1Offset;
		int s2Offset;
		if (start2RelativeToStart1 < 0) {
			s1Offset = 0;
			s2Offset = -start2RelativeToStart1;
		} else {
			s1Offset = start2RelativeToStart1;
			s2Offset = 0;
		}
		for (int i = 0; i + s1Offset < s1.baseCount && i + s2Offset < s2.baseCount; i += BASES_PER_WORD) {
			long s1bases = s1.get2BitBases(i + s1Offset);
			long s2bases = s2.get2BitBases(i + s2Offset);
			// match pairs of bits
			long matchVector = ~(s1bases ^ s2bases);
			matchVector &= matchVector >>> 1;
			matchVector &= MATCH_MASK; // only keep the matches within each base (/bit pair)
			// zero out LSBs running off the end of either sequence
			int s1remainingLength = s1.baseCount - i - s1Offset;
			int s2remainingLength = s2.baseCount - i - s2Offset;
			int basesToConsider = Math.min(BASES_PER_WORD, Math.min(s1remainingLength, s2remainingLength));
			matchVector &= -1L << ((BASES_PER_WORD - basesToConsider) * BITS_PER_BASE);
			matches += Long.bitCount(matchVector);
		}
		return matches;
	}
	public static int overlapLength(PackedSequence s1, PackedSequence s2, int start2RelativeToStart1) {
		return IntervalUtil.overlapsWidthClosed(0, s1.baseCount - 1, start2RelativeToStart1, start2RelativeToStart1 + s2.baseCount - 1);
	}

	/**
	 *
	 * @param psFixed sequence to compare starting at the first base
	 * @param psOffset sequence to compare starting at the offsetBases base
	 * @param offsetBases number of bases to offset the comparison from the start of the fixed sequence. Cannot be negative
	 * @return
	 */
	private static int optimised_overlapMatches(PackedSequence psFixed, PackedSequence psOffset, int offsetBases) {
		// TODO: Java support for SIMD
		int overlapLength = psFixed.baseCount - offsetBases;
		assert (offsetBases >= 0);
		assert (overlapLength <= 0);
		if (overlapLength <= BASES_PER_WORD) {
			// Less than a single word to compare
			long bases = psOffset.get2BitBases(offsetBases);
			long matchVector = ~(psFixed.packed[0] ^ bases);
			matchVector &= matchVector >>> 1;
			long mask = MATCH_MASK;
			mask <<= 2 * (BASES_PER_WORD - overlapLength); // ignore the LSB bases
			matchVector &= mask;
			return Long.bitCount(matchVector);
		}
		while (overlapLength >= BASES_PER_WORD) {
			// TODO full word comparison
		}
		// TODO compare remaining bases in final word
		throw new IllegalStateException("NYI");
	}
	// TODO: optimise:
	// - function is symmetric so arg order doesn't matter
	// - fix one sequence, loop over the full length
	//
	/**
	 * Returns the 32 bases starting with the base at the given offset
	 * @param start 0-based offset of sequence to extract
	 * @return packed sequenced in MSBs
	 */
	private long get2BitBases(int start) {
		assert(start < baseCount);
		assert(start >= 0);
		int wordIndex = start / BASES_PER_WORD;
		int wordBaseIndex = start % BASES_PER_WORD;
		long result = packed[wordIndex] << (wordBaseIndex * BITS_PER_BASE);
		if (wordIndex + 1 < packed.length && wordBaseIndex > 0) {
			// grab second word
			result |= packed[wordIndex + 1] >>> ((BASES_PER_WORD - wordBaseIndex) * BITS_PER_BASE);
		}
		// TODO: zero out bases after length?
		return result;
	}
	public int length() {
		return baseCount;
	}
	public long[] asLongArray() {
		return packed;
	}

    public void setKmer(long kmer, int offset, int k) {
		long existingKmer = getKmer(offset, k);
		for (int i = 0; i < k; i++) {
			long existingBase = (existingKmer >> (BITS_PER_BASE * i)) & BASE_MASK;
			long newBase = (kmer >> (BITS_PER_BASE * i)) & BASE_MASK;
			if (existingBase != newBase) {
				// first base is in the MSBs of the kmer
				setBaseEncoded(offset + k - 1 - i, newBase);
			}
		}
    }
}