package au.edu.wehi.idsv.visualisation;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.debruijn.subgraph.SubgraphPathNode;

/**
 * Tracks de Bruijn subgraph assembly progress for a single subgraph
 * @author cameron.d
 *
 */
public interface SubgraphAssemblyAlgorithmTracker {

	public abstract void finalAnchors(int minAnchor, int maxAnchor);

	public abstract void splitOutReferencePaths(int pathsSplits);

	public abstract void calcNonReferenceSubgraphs(List<Set<SubgraphPathNode>> subgraphs, List<Set<SubgraphPathNode>> startingNodes);

	public abstract void assemblyNonReferenceContigs(List<List<SubgraphPathNode>> assembledPaths, List<LinkedList<Long>> assembledKmers, int nodesTraversed);

	public abstract void toAssemblyEvidence(VariantContextDirectedEvidence assembly);

	public abstract void shrink(int edgesRemoved, int nodesCollapsed);

	public abstract void collapse(int collapseIterations, int pathsCollapsed, int leavesCollapsed, int nodesTraversed, int nodeCountReducedBy);

	public abstract void assemblyStarted();

	public abstract void assemblyComplete();

	public abstract void generatePathGraph(int kmers, int nodes, int edges);

	public abstract String toBed();
}