package au.edu.wehi.idsv.debruijn.positional;

import java.util.List;

public class PathGraphNode {
	private int startPosition;
	private int endPosition;
	private int weight;
	private List<Long> kmers;
}
