package au.edu.wehi.idsv.debruijn;

public class ReadKmer {
	public ReadKmer(long kmer, long weight) {
		this.kmer = kmer;
		this.weight = weight;
	}
	public final long kmer;
	public final long weight;
	/**
	 * Display as much of the kmer as we can
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		long state = kmer;
		while (state > 0) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(state, 0));
			state >>=2;
		}
		return sb.reverse().toString();
	}
	public String toString(int k) {
		StringBuilder sb = new StringBuilder();
		long state = kmer;
		for (int i = 0; i < k; i++) {
			sb.append((char)KmerEncodingHelper.lastBaseEncodedToPicardBase(state, 0));
			state >>=2;
		}
		return sb.reverse().toString();
	}
}
