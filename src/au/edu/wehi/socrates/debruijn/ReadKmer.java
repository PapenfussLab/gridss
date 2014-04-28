package au.edu.wehi.socrates.debruijn;

public class ReadKmer {
	public ReadKmer(long kmer, long weight) {
		this.kmer = kmer;
		this.weight = weight;
	}
	public final long kmer;
	public final long weight;
}
