package au.edu.wehi.idsv.debruijn;

import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import com.google.common.primitives.Bytes;


public class KmerEncodingHelper {
	private KmerEncodingHelper() { }
	/**
	 * Maximum kmer size able to be encoded in a long
	 */
	public static final int MAX_K = Long.SIZE / 2;
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
	private static final long[] complementBits = new long[MAX_K + 1];
	/**
	 * Bit mask return only bits used by a kmer of length k 
	 */
	private static final long[] usedBits = new long[MAX_K + 1];
	private static final byte[] ENCODED_TO_PICARD_LOOKUP = { 'T', 'C', 'A', 'G' };
	static {
		long usedMask = 0;
		long complementMask = 0;
		for (int i = 0; i <= MAX_K; i++) {
			complementBits[i] = complementMask;
			complementMask <<= 2;
			complementMask |= 2;
			usedBits[i] = usedMask;
			usedMask <<= 2;
			usedMask |= 3;
		}
	}
	/**
	 * Converts a base to 2bit representation
	 * @see http://genome.ucsc.edu/FAQ/FAQformat#format7
	 * @param base ASCII byte representation of base
	 * @return 2bit representation of base
	 */
	public static int picardBaseToEncoded(byte base) {
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
	public static boolean isAmbiguous(byte base) {
		switch (base) {
			case 'G':
			case 'g':
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'T':
			case 't':
				return false;
			default:
				return true;
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
	public static byte encodedToPicardBase(int encoded) {
		return ENCODED_TO_PICARD_LOOKUP[(int)encoded & 0x03];
		//switch ((int)encoded & 0x03) {
		//	case 3: return 'G';
		//	case 1: return 'C';
		//	case 0: return 'T';
		//	default:
		//	case 2: return 'A';
		//}
	}
	public static byte encodedToPicardBase(long encoded) {
		return encodedToPicardBase((int)encoded);
	}
	public static byte firstBaseEncodedToPicardBase(int k, long state) {
		return encodedToPicardBase(state >>> (2 * (k - 1)));
	}
	public static byte lastBaseEncodedToPicardBase(long state) {
		return encodedToPicardBase(state);
	}
	public static void assertValid(int k, long encoded) {
		if ((encoded & ~usedBits[k]) != 0) {
			throw new IllegalArgumentException(String.format("Sanity check failure: state %d is not a %dmer", encoded, k));
		}
	}
	/**
	 * Converts a 2bit encoded kmer into picard bases  
	 * @param base
	 * @see http://genome.ucsc.edu/FAQ/FAQformat#format7
	 * @return 2bit
	 */
	public static byte[] encodedToPicardBases(int k, long encoded) {
		assertValid(k, encoded);
		long state = encoded;
		byte[] result = new byte[k];
		for (int i = 0; i < k; i++) {
			result[k - i - 1] = lastBaseEncodedToPicardBase(state);
			state >>>= 2;
		}
		return result;
	}
	public static boolean isPrev(int k, long encoded, long prev) {
		return isNext(k, prev, encoded);
	}
	public static boolean isNext(int k, long encoded, long next) {
		return clearBase(k - 1, encoded) == next >>> 2;
	}
	public static long[] nextStates(int k, long encoded) {
		long next = clearBase(k - 1, encoded) << 2;
		return new long[] {
				next,
				next | 1L,
				next | 2L,
				next | 3L
		};
	}
	public static long[] prevStates(int k, long encoded) {
		long next = encoded >>> 2;
		return new long[] {
				next,
				next | (1L << (2 * k - 2)),
				next | (2L << (2 * k - 2)),
				next | (3L << (2 * k - 2))
		};
	}
	/**
	 * returns all kmers adjacent to this kmer
	 * @param k k
	 * @param encoded kmer
	 * @return subsequent kmers, followed by previous kmers
	 */
	public static long[] adjacentStates(int k, long encoded) {
		long next = clearBase(k - 1, encoded) << 2;
		long prev = encoded >>> 2;
		return new long[] {
				next,
				next | 1L,
				next | 2L,
				next | 3L,
				prev,
				prev | (1L << (2 * k - 2)),
				prev | (2L << (2 * k - 2)),
				prev | (3L << (2 * k - 2))
		};
	}
	public static long nextState(int k, long state, byte picardBase) {
		assertValid(k, state);
		long next = clearBase(k - 1, state) << 2;
		next = next | picardBaseToEncoded(picardBase);
		assertValid(k, next);
		return next;
	}
	public static boolean lastBaseMatches(int k, long state1, long state2) {
		return (state1 & 3L) == (state2 & 3L);
	}
	public static boolean firstBaseMatches(int k, long state1, long state2) {
		return state1 >>> ((k - 1) * 2) == state2 >>> ((k - 1) * 2);
	}
	private static long clearBase(int k, long state) {
		long bitsToClear = (1L << ((2*k)+1)) | (1L << (2*k));
		return state & ~bitsToClear;
	}
	public static String toString(int k, long state) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < k; i++) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(state));
			state >>>=2;
		}
		return sb.reverse().toString();
	}
	/**
	 * Converts the state to the best string with unknown k.
	 * @param state encoded kmer
	 * @return base sequence missing leading Ts.
	 */
	public static String toApproximateString(long state) {
		int nonZeroPos = 64 - Long.numberOfLeadingZeros(state);
		return toString(Math.max((nonZeroPos + 1) / 2, 1), state);
	}
	/**
	 * Returns the number of bases difference between the two states
	 * @param k
	 * @param state1 kmer
	 * @param state2 kmer
	 * @return number of mismatching bases
	 */
	public static int basesDifference(int k, long state1, long state2) {
		return k - basesMatching(k, state1, state2);
	}
	/**
	 * Returns the number of bases matching between the two states
	 * @param state1 kmer
	 * @param state2 kmer
	 * @return number of bases identical in both kmers 
	 */
	public static int basesMatching(int k, long state1, long state2) {
		// set lower bit to 1 if both high and low bits match
		long bitsMatch = ~(state1 ^ state2);
		long baseMatch = ((bitsMatch & HIGH_BITS) >>> 1) & bitsMatch;
		int baseCount = Long.bitCount(baseMatch & usedBits[k]); 
		return baseCount;
	}
	/**
	 * Base calls of contig
	 * @param path kmer contig
	 * @return base calls of a positive strand SAMRecord readout of contig
	 */
	public static byte[] baseCalls(List<Long> path, int k) {
		int assemblyLength = path.size() + k - 1;
		byte[] bases = KmerEncodingHelper.encodedToPicardBases(k, path.get(0));
		bases = Arrays.copyOf(bases, assemblyLength);
		int offset = k - 1;
		for (Long node : path) {
			bases[offset] = KmerEncodingHelper.lastBaseEncodedToPicardBase(node);
			offset++;
		}
		return bases;
	}
	/**
	 * Base calls of contig
	 * @param path kmer contig
	 * @return base calls of a positive strand SAMRecord readout of contig
	 */
	public static int totalBaseDifference(Iterator<Long> pathA, Iterator<Long> pathB, int k) {
		if (pathA == null || pathB == null) return 0;
		if (!pathA.hasNext() || !pathB.hasNext()) return 0;
		int diffCount = basesDifference(k, pathA.next(), pathB.next());
		while (pathA.hasNext() && pathB.hasNext()) {
			if (!lastBaseMatches(k, pathA.next(), pathB.next())) {
				diffCount++;
			}
		}
		return diffCount;
	}
	/**
	 * Sums base counts for the given sequence
	 * @return
	 */
	public static int[] baseCounts(int k, LongArrayList path) {
		int[] counts = new int[4];
		long startKmer = path.getLong(0);
		for (int i = 0; i < k; i++) {
			counts[(int)startKmer & 3]++;
			startKmer >>>= 2;
		}
		for (int i = 1; i < path.size(); i++) {
			counts[(int)path.getLong(i) & 3]++;
		}
		return counts;
	}
}
