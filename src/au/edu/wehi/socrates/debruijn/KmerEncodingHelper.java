package au.edu.wehi.socrates.debruijn;

import java.math.BigInteger;

public class KmerEncodingHelper {
	private KmerEncodingHelper() { }
	/**
	 * Converts a base to 2bit representation
	 * @see http://genome.ucsc.edu/FAQ/FAQformat#format7
	 * @param base
	 * @return 2bit
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
	public static long picardBaseToEncoded(int k, byte[] bases) {
		if (bases == null) throw new NullPointerException("bases null");
		if (k > bases.length) throw new IllegalArgumentException("fewer bases than k");
		long result = picardBaseToEncoded(bases[0]);
		for (int i = 1; i < k; i++) {
			result <<= 2;
			result |= picardBaseToEncoded(bases[i]);
		}
		return result;
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
		return encodedToPicardBase(state >> (2 * (k - 1)));
	}
	public static byte lastBaseEncodedToPicardBase(long state, int k) {
		return encodedToPicardBase(state);
	}
	/**
	 * Converts a 2bit encoded kmer into picard bases  
	 * @see http://genome.ucsc.edu/FAQ/FAQformat#format7
	 * @param base
	 * @return 2bit
	 */
	public static byte[] encodedToPicardBases(long encoded, int length) {
		byte[] result = new byte[length];
		for (int i = 0; i < length; i++) {
			result[length - i - 1] = lastBaseEncodedToPicardBase(encoded, length);
			encoded >>= 2;
		}
		return result;
	}
	public static long[] nextStates(long encoded, int length) {
		long next = clearBase(encoded, length - 1) << 2;
		return new long[] {
				next,
				next | 1,
				next | 2,
				next | 3
		};
	}
	public static long[] prevStates(long encoded, int length) {
		long next = encoded >> 2;
		return new long[] {
				next,
				next | (1 << (2 * length - 2)),
				next | (2 << (2 * length - 2)),
				next | (3 << (2 * length - 2))
		};
	}
	public static long nextState(int k, long state, byte picardBase) {
		long next = clearBase(state, k - 1) << 2;
		return next | picardBase;
	}
	private static long clearBase(long state, int offset) {
		return state & ~((1 << (2*offset+1)) & (1 << (2*offset)));
	}
}
