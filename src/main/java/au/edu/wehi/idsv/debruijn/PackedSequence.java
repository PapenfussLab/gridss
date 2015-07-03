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
	private static final int BASES_PER_WORD = Integer.SIZE / BITS_PER_BASE;
	/**
	 * First base is packed in MSB of first word
	 * Second base is packed in second MSB of first word
	 * and so on
	 */
	protected final int[] packed;
	protected static int[] pack(byte[] bases, boolean reverse, boolean complement) {
		int[] packed = new int[IntMath.divide(bases.length, BASES_PER_WORD, RoundingMode.CEILING)];
		int currentWord = 0;
		int packedIntoCurrentWord = 0;
		int currentOffset = 0;
		for (int i = 0; i < bases.length; i++) {
			byte base = bases[reverse ? bases.length - 1 - i: i];
			long encoded = KmerEncodingHelper.picardBaseToEncoded(base);
			if (complement) {
				encoded = KmerEncodingHelper.complement(1, encoded);
			}
			currentWord <<= 2;
			currentWord |= (int)encoded;
			packedIntoCurrentWord ++;
			if (packedIntoCurrentWord == BASES_PER_WORD) {
				packed[currentOffset++] = currentWord;
				packedIntoCurrentWord = 0;
				currentWord = 0;
			}
		}
		if (packedIntoCurrentWord != 0) {
			packed[currentOffset] = (int)(currentWord << (BITS_PER_BASE * (BASES_PER_WORD - packedIntoCurrentWord)));
		}
		return packed;
	}
	public PackedSequence(byte[] bases, boolean reverse, boolean complement) {
		this.packed = pack(bases, reverse, complement);
	}
	public byte get(int offset) {
		int wordIndex = (offset) / BASES_PER_WORD;
		int wordOffset = BASES_PER_WORD - ((offset) % BASES_PER_WORD) - 1;
		byte b = KmerEncodingHelper.encodedToPicardBase(packed[wordIndex] >>> (wordOffset * BITS_PER_BASE));
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
	public long getEncodedLong(int offset, int length) {
		assert(length <= Long.SIZE / BITS_PER_BASE);
		assert(offset + length <= packed.length * BASES_PER_WORD);
		int toTake = length;
		long accum = 0;
		while (toTake > 0) {
			int wordOffset = offset / BASES_PER_WORD;
			int basesToSkipInWord = offset % BASES_PER_WORD;
			int basesToTakeInWord = Math.min(toTake, BASES_PER_WORD - basesToSkipInWord);
			
			accum <<= basesToTakeInWord * BITS_PER_BASE;
			accum |= ((int)packed[wordOffset] >> (BITS_PER_BASE * (BASES_PER_WORD - basesToSkipInWord - basesToTakeInWord))) & ((1 << (BITS_PER_BASE * basesToTakeInWord)) - 1);
			toTake -= basesToTakeInWord;
			offset += basesToTakeInWord;
		}
		return accum;
	}
}