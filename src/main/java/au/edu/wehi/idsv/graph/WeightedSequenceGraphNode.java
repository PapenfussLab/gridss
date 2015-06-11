package au.edu.wehi.idsv.graph;

import java.util.List;

public interface WeightedSequenceGraphNode {
	int length();
	int weight();
	int weight(int offset);
}
