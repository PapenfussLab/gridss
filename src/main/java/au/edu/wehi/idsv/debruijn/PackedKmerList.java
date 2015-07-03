package au.edu.wehi.idsv.debruijn;


public class PackedKmerList extends PackedSequence {
	private final byte[] weights;
	private final byte k;
	public PackedKmerList(int k, byte[] bases, byte[] qual, boolean reverse, boolean complement) {
		super(bases, reverse, complement);
		int kmers = bases.length - k + 1;
		this.k = (byte)k;
		if (kmers <= 0) {
			this.weights = new byte[0];
		} else if (qual == null) {
			this.weights = new byte[kmers];
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
		return getEncodedLong(offset, k); 
	}
	public int weight(int offset) {
		return weights[offset];
	}
	public int length() {
		return weights.length;
	}
}
