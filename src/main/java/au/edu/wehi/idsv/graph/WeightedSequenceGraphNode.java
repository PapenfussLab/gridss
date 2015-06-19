package au.edu.wehi.idsv.graph;


public interface WeightedSequenceGraphNode {
	int length();
	int weight();
	int weight(int offset);
}
