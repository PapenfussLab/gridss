package au.edu.wehi.idsv.debruijn;

import java.math.RoundingMode;

import com.google.common.math.IntMath;

/**
 * Compresses the given sequence by representing in 2-bit format
 * @author Daniel Cameron
 *
 */
public class PackedSequence {
	private static final int BITS_PER_BASE = 2;
	private static final int BASES_PER_WORD = Long.SIZE / BITS_PER_BASE;
	//private static final long BASE_MASK = (1 << BITS_PER_BASE) - 1;
	private static final int ARRAY_SHIFT = Long.SIZE - 1 - Long.numberOfLeadingZeros(BASES_PER_WORD);
	private static final int ARRAY_OFFSET_MASK = (1 << ARRAY_SHIFT) - 1;
	private static final long BASE_MASK = (1 << BITS_PER_BASE) - 1;
	/**
	 * First base is packed in MSB of first word
	 * Second base is packed in second MSB of first word
	 * and so on
	 */
	private final long[] packed;
	public PackedSequence(byte[] bases, boolean reverse, boolean complement) {
		packed = new long[IntMath.divide(bases.length, BASES_PER_WORD, RoundingMode.CEILING)];
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
	private void setBaseEncoded(final int offset, final long base) {
		int wordIndex = offset >> ARRAY_SHIFT;
		int wordOffset = BASES_PER_WORD - 1 - (offset & ARRAY_OFFSET_MASK);
		assert(wordIndex < packed.length);
		long word = packed[wordIndex];
		word &= ~(BASE_MASK << (BITS_PER_BASE * wordOffset));
		word |= base << (BITS_PER_BASE * wordOffset);
		packed[wordIndex] = word;
	}
	private long getBaseEncoded(final int offset) {
		int wordIndex = offset >> ARRAY_SHIFT;
		int wordOffset = BASES_PER_WORD - 1 - (offset & ARRAY_OFFSET_MASK);
		assert(wordIndex < packed.length);
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
		return new String(getBytes(0, packed.length * BASES_PER_WORD));
	}
}