package au.edu.wehi.idsv.visualisation;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.debruijn.subgraph.SubgraphPathNode;

/**
 * Does not perform any algorithm tracking
 * @author cameron.d
 *
 */
public class NontrackingSubgraphTracker implements SubgraphAssemblyAlgorithmTracker {
	@Override
	public void finalAnchors(int minAnchor, int maxAnchor) {
	}
	@Override
	public void splitOutReferencePaths(int pathsSplits) {
	}
	@Override
	public void calcNonReferenceSubgraphs(
			List<Set<SubgraphPathNode>> subgraphs,
			List<Set<SubgraphPathNode>> startingNodes) {
	}
	@Override
	public void assemblyNonReferenceContigs(
			List<List<SubgraphPathNode>> assembledPaths,
			List<LinkedList<Long>> assembledKmers, int nodesTraversed) {
	}
	@Override
	public void toAssemblyEvidence(VariantContextDirectedEvidence assembly) {
	}
	@Override
	public void shrink(int edgesRemoved, int nodesCollapsed) {
	}
	@Override
	public void assemblyStarted() {
	}
	@Override
	public void assemblyComplete() {
	}
	@Override
	public void generatePathGraph(int nodes, int edges) {
	}
	@Override
	public void collapse(int collapseIterations, int pathsCollapsed,
			int leavesCollapsed, int nodesTraversed, int nodeCountReducedBy) {
		// TODO Auto-generated method stub
		
	}
}
