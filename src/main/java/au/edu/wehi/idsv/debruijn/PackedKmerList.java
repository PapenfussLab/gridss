package au.edu.wehi.idsv.debruijn;


public class PackedKmerList {
	private final byte[] weights;
	protected final byte k;
	protected PackedSequence seq;
	public PackedKmerList(int k, byte[] bases, byte[] qual, boolean reverse, boolean complement) {
		this.seq = new PackedSequence(bases, reverse, complement);
		this.k = (byte)k;
		if (length() <= 0 || qual == null) {
			this.weights = null;
		} else {
			this.weights = calcWeight(k, qual, reverse);
		}
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
		return seq.getKmer(offset, k);
	}
	public int weight(int offset) {
		assert(offset < length());
		if (weights == null) return 0;
		return weights[offset];
	}
	public int length() {
		return seq.kmers(k);
	}
	public int kmerSize() {
		return k;
	}
}
