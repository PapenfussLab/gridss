package au.edu.wehi.idsv.debruijn;

import java.util.List;

import com.google.common.primitives.Bytes;


public class KmerEncodingHelper {
	private KmerEncodingHelper() { }
	/**
	 * Maximum kmer size able to be encoded in a long
	 */
	public static final int MAX_KMER = Long.SIZE / 2;
	/**
	 * Every high bit of each 2bit base
	 */
	private static final long HIGH_BITS = 0xAAAAAAAAAAAAAAAAL;
	/**
	 * Every low bit of each 2bit base
	 */
	private static final long LOW_BITS = 0x5555555555555555L;
	/**
	 * Bits patterns required to complement an encoded kmer of length k.
	 * Conveniently, complementing a UCSC 2bit encoded base requires just flipping the high bit.
	 */
	private static final long[] complementBits = new long[MAX_KMER + 1];
	static {
		long pattern = 0;
		for (int i = 0; i <= MAX_KMER; i++) {
			complementBits[i] = pattern;
			pattern <<= 2;
			pattern |= 2;
		}
	}
	/**
	 * Converts a base to 2bit representation
	 * @see http://genome.ucsc.edu/FAQ/FAQformat#format7
	 * @param base ASCII byte representation of base
	 * @return 2bit representation of base
	 */
	private static long picardBaseToEncoded(byte base) {
		switch (base) {
			case 'G':
			case 'g':
				return 3;
			case 'A':
			case 'a':
				return 2;
			case 'C':
			case 'c':
				return 1;
			case 'T':
			case 't':
				return 0;
			default:
				// default to A for unknown bases
				return 2;
		}
	}
	public static long picardBaseToEncoded(int k, List<Byte> bases) {
		if (bases == null) throw new NullPointerException("bases null");
		if (k > bases.size()) throw new IllegalArgumentException("fewer bases than k");
		long result = picardBaseToEncoded(bases.get(0));
		for (int i = 1; i < k; i++) {
			result <<= 2;
			result |= picardBaseToEncoded(bases.get(i));
		}
		assertValid(k, result);
		return result;
	}
	public static long picardBaseToEncoded(int k, byte[] bases) {
		return picardBaseToEncoded(k, Bytes.asList(bases));
	}
	/**
	 * Reverses the kmer bases
	 * @param k kmer size
	 * @param encoded encoded kmer
	 * @return reversed kmer
	 */
	public static long reverse(int k, long encoded) {
		// independently reverse each of the high and low bits of each base
		long reversed = (Long.reverse(encoded & HIGH_BITS) << 1) | (Long.reverse(encoded & LOW_BITS) >>> 1);
		// and shuffle back so the bits we're using are the lowest 2k bits
		return reversed >>> (Long.SIZE - 2*k);
	}
	public static long complement(int k, long encoded) {
		return encoded ^ complementBits[k];
	}
	public static long reverseComplement(int k, long encoded) {
		return reverse(k, complement(k, encoded));
	}
	private static byte encodedToPicardBase(long encoded) {
		switch ((int)encoded & 0x03) {
			case 3: return 'G';
			case 1: return 'C';
			case 0: return 'T';
			default:
			case 2: return 'A';
		}
	}
	public static byte firstBaseEncodedToPicardBase(long state, int k) {
		return encodedToPicardBase(state >>> (2 * (k - 1)));
	}
	public static byte lastBaseEncodedToPicardBase(long state, int k) {
		return encodedToPicardBase(state);
	}
	public static void assertValid(int k, long encoded) {
		if (k != 32 && encoded >>> (2 * k) != 0) {
			throw new IllegalArgumentException(String.format("Sanity check failure: state %d is not a %dmer", encoded, k));
		}
	}
	/**
	 * Converts a 2bit encoded kmer into picard bases  
	 * @see http://genome.ucsc.edu/FAQ/FAQformat#format7
	 * @param base
	 * @return 2bit
	 */
	public static byte[] encodedToPicardBases(long encoded, int length) {
		assertValid(length, encoded);
		long state = encoded;
		byte[] result = new byte[length];
		for (int i = 0; i < length; i++) {
			result[length - i - 1] = lastBaseEncodedToPicardBase(state, length);
			state >>>= 2;
		}
		return result;
	}
	public static long[] nextStates(long encoded, int length) {
		long next = clearBase(length - 1, encoded) << 2;
		return new long[] {
				next,
				next | 1,
				next | 2,
				next | 3
		};
	}
	public static long[] prevStates(long encoded, int length) {
		long next = encoded >>> 2;
		return new long[] {
				next,
				next | (1 << (2 * length - 2)),
				next | (2 << (2 * length - 2)),
				next | (3 << (2 * length - 2))
		};
	}
	public static long nextState(int k, long state, byte picardBase) {
		assertValid(k, state);
		long next = clearBase(k - 1, state) << 2;
		next = next | picardBaseToEncoded(picardBase);
		assertValid(k, next);
		return next;
	}
	private static long clearBase(int k, long state) {
		long bitsToClear = (1L << ((2*k)+1)) | (1L << (2*k));
		return state & ~bitsToClear;
	}
	public static String toString(int k, long state) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < k; i++) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(state, 0));
			state >>>=2;
		}
		return sb.reverse().toString();
	}
}
