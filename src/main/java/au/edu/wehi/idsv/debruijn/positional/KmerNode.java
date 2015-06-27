package au.edu.wehi.idsv.debruijn.positional;


public interface KmerNode {
	long kmer();
	int startPosition();
	int endPosition();
	int weight();
	boolean isReference();
	default int width() { return endPosition() - startPosition() + 1; }
	default int firstKmerStartPosition() { return startPosition(); }
	default int firstKmerEndPosition() { return endPosition(); }
	default long firstKmer() { return kmer(); }
	default int length() { return 1; }
}
