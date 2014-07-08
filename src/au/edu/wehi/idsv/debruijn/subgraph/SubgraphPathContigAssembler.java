package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.debruijn.PathNode;

/**
 * Reduces the de bruijn graph to a set of paths before assembling
 * @author Daniel Cameron
 *
 */
public abstract class SubgraphPathContigAssembler implements SubgraphContigAssembler {
	protected final DeBruijnReadGraph graph;
	protected final DeBruijnPathGraph<DeBruijnSubgraphNode, SubgraphPathNode> pathGraph;
	protected final long seedKmer;
	public SubgraphPathContigAssembler(DeBruijnReadGraph graph, long seedKmer) {
		this.graph = graph;
		this.seedKmer = seedKmer;
		this.pathGraph = new SubgraphPathGraph(graph, seedKmer);
	}
	@Override
	public abstract List<LinkedList<Long>> assembleContigs();
	//pathGraph.removeSelfIntersectingPaths();
	//pathGraph.collapseSimilarPaths(2, true);
	// return only the best path on the subgraph
}
