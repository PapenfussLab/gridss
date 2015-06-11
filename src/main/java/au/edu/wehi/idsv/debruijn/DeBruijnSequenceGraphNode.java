package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.graph.WeightedSequenceGraphNode;


public interface DeBruijnSequenceGraphNode extends WeightedSequenceGraphNode {
	long kmer(int offset);
}
