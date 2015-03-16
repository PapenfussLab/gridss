package au.edu.wehi.idsv.visualisation;

import java.util.LinkedList;
import java.util.List;

import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.debruijn.subgraph.SubgraphPathNode;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

/**
 * Tracks de Bruijn subgraph assembly progress for a single subgraph
 * @author cameron.d
 *
 */
public interface SubgraphAssemblyAlgorithmTracker {

	public abstract void finalAnchors(int minAnchor, int maxAnchor);

	public abstract void splitOutReferencePaths(int pathsSplits);

	//public abstract void calcNonReferenceSubgraphs(List<Set<SubgraphPathNode>> subgraphs, List<Set<SubgraphPathNode>> startingNodes);

	public abstract void assemblyNonReferenceContigs(List<List<SubgraphPathNode>> assembledPaths, List<LinkedList<Long>> assembledKmers, int nodesTraversed);

	public abstract void toAssemblyEvidence(VariantContextDirectedEvidence assembly);

	public abstract void shrink(int edgesRemoved, int nodesCollapsed);

	public abstract void collapse(int collapseIterations, int pathsCollapsed, int leavesCollapsed, int nodesTraversed, int nodeCountReducedBy);

	public abstract void assemblyStarted();

	public abstract void assemblyComplete();

	public abstract void generatePathGraph(int kmers, int nodes, int edges);

	public abstract String toBed();
	
	public abstract int getReferenceIndex();
	
	public abstract long getStartAnchorPosition();
	
	public abstract long getEndAnchorPosition();
	
	public static Ordering<SubgraphAssemblyAlgorithmTracker> ByGenomicPosition = new Ordering<SubgraphAssemblyAlgorithmTracker>() {
		@Override
		public int compare(SubgraphAssemblyAlgorithmTracker left, SubgraphAssemblyAlgorithmTracker right) {
			return ComparisonChain.start()
			        .compare(left.getReferenceIndex(), right.getReferenceIndex())
			        .compare(left.getStartAnchorPosition(), right.getStartAnchorPosition())
			        .compare(left.getEndAnchorPosition(), right.getEndAnchorPosition())
			        .result();
		}
	};
}