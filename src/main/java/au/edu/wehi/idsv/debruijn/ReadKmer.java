package au.edu.wehi.idsv.debruijn;

public class ReadKmer {
	public ReadKmer(long kmer, int weight, boolean containsAmbiguousBases) {
		this.kmer = kmer;
		this.weight = weight;
		this.containsAmbiguousBases = containsAmbiguousBases;
	}
	public final long kmer;
	public final int weight;
	public final boolean containsAmbiguousBases;
	/**
	 * Display as much of the kmer as we can
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		long state = kmer;
		while (state > 0) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(state));
			state >>=2;
		}
		return sb.reverse().toString();
	}
	public String toString(int k) {
		StringBuilder sb = new StringBuilder();
		long state = kmer;
		for (int i = 0; i < k; i++) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(state));
			state >>=2;
		}
		return sb.reverse().toString();
	}
}
