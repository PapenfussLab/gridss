package au.edu.wehi.idsv.visualisation;

import java.util.List;

import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.graph.PathNode;

/**
 * Does not perform any algorithm tracking
 * @author Daniel Cameron
 *
 */
public class NontrackingSubgraphTracker<T, PN extends PathNode<T>> implements SubgraphAssemblyAlgorithmTracker<T, PN> {
	@Override
	public void finalAnchors(int minAnchor, int maxAnchor) {
	}
	@Override
	public void splitOutReferencePaths(int pathsSplits) {
	}
	@Override
	public void assemblyNonReferenceContigs(List<List<PN>> assembledPaths, int nodesTraversed) {
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
	public void collapse(int collapseIterations, int pathsCollapsed,
			int leavesCollapsed, int nodesTraversed, int nodeCountReducedBy) {
	}
	@Override
	public String toBed() {
		return null;
	}
	@Override
	public void generatePathGraph(int kmers, int nodes, int edges) {
	}
	@Override
	public int getReferenceIndex() {
		return -1;
	}
	@Override
	public long getStartAnchorPosition() {
		return 0;
	}
	@Override
	public long getEndAnchorPosition() {
		return 0;
	}
}
