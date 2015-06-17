package au.edu.wehi.idsv.debruijn;

import java.math.RoundingMode;

import com.google.common.math.IntMath;


public class PackedKmerList {
	private static final int BITS_PER_BASE = 2;
	private static final int BASES_PER_WORD = Byte.SIZE / BITS_PER_BASE;
	private final byte[] twoBitBases;
	private final byte[] weights;
	private final byte k;
	public PackedKmerList(int k, byte[] bases, byte[] qual, boolean reverse, boolean complement) {
		int kmers = bases.length - k + 1;
		this.k = (byte)k;
		this.twoBitBases = twoBitEncodeBases(k, bases, reverse, complement);
		if (kmers <= 0) {
			this.weights = new byte[0];
		} else if (qual == null) {
			this.weights = new byte[kmers];
		} else {
			this.weights = calcWeight(k, qual, reverse);
		}
	}
	private static byte[] twoBitEncodeBases(int k, byte[] bases, boolean reverse, boolean complement) {
		byte[] kmer = new byte[IntMath.divide(bases.length, Byte.SIZE / 2, RoundingMode.CEILING)];
		byte currentWord = 0;
		int packedIntoCurrentWord = 0;
		int currentOffset = 0;
		for (int i = 0; i < bases.length; i++) {
			byte base = bases[reverse ? bases.length - 1 - i: i];
			long encoded = KmerEncodingHelper.picardBaseToEncoded(base);
			if (complement) {
				encoded = KmerEncodingHelper.complement(1, encoded);
			}
			currentWord <<= 2;
			currentWord |= (byte)encoded;
			packedIntoCurrentWord ++;
			if (packedIntoCurrentWord == BASES_PER_WORD) {
				kmer[currentOffset++] = currentWord;
				packedIntoCurrentWord = 0;
				currentWord = 0;
			}
		}
		if (packedIntoCurrentWord != 0) {
			kmer[currentOffset] = (byte)(currentWord << (BITS_PER_BASE * (BASES_PER_WORD - packedIntoCurrentWord)));
		}
		return kmer;
	}
	private static byte[] calcWeight(int k, byte[] qual, boolean reverse) {
		byte[] weights = new byte[qual.length - k + 1];
		// TODO: see if keeping a small lookup tree is faster
		for (int i = 0; i < weights.length; i++) {
			weights[i] = Byte.MAX_VALUE;
		}
		for (int i = 0; i < qual.length; i++) {
			// +1 to ensure all weights are non-zero
			int weight = qual[reverse ? qual.length - 1 - i : i] + 1;
			for (int j = -(k - 1); j <= 0; j++) {
				int pos = i + j;
				if (pos >= 0 && pos < weights.length) {
					weights[pos] = (byte)Math.min(weights[pos], weight);
				}
			}
		}
		return weights;
	}
	public long kmer(int offset) {
		int toTake = k;
		long accum = 0;
		while (toTake > 0) {
			int wordOffset = offset / BASES_PER_WORD;
			int basesToSkipInWord = offset % BASES_PER_WORD;
			int basesToTakeInWord = Math.min(toTake, BASES_PER_WORD - basesToSkipInWord);
			
			accum <<= basesToTakeInWord * BITS_PER_BASE;
			accum |= ((int)twoBitBases[wordOffset] >> (BITS_PER_BASE * (BASES_PER_WORD - basesToSkipInWord - basesToTakeInWord))) & ((1 << (BITS_PER_BASE * basesToTakeInWord)) - 1);
			toTake -= basesToTakeInWord;
			offset += basesToTakeInWord;
		}
		return accum;
	}
	public int weight(int offset) {
		return weights[offset];
	}
	public int length() {
		return weights.length;
	}
}
