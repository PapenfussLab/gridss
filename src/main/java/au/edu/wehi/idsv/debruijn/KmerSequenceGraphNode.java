package au.edu.wehi.idsv.debruijn;

public interface KmerSequenceGraphNode {
	int length();
	int weight();
	int weight(int offset);
}
