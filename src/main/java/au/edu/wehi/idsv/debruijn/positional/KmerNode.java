package au.edu.wehi.idsv.debruijn.positional;


public interface KmerNode {
	long lastKmer();
	int lastStart();
	int lastEnd();
	int weight();
	boolean isReference();
	default int width() { return lastEnd() - lastStart() + 1; }
	default int firstKmerStart() { return lastStart(); }
	default int firstKmerEnd() { return lastEnd(); }
	default long firstKmer() { return lastKmer(); }
	default int length() { return 1; }
}
