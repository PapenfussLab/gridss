package au.edu.wehi.idsv.visualisation;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.debruijn.subgraph.SubgraphPathNode;

import com.google.common.collect.Lists;

/**
 * Tracks algorithm metrics for each de Bruijn subgraph assembly 
 * @author cameron.d
 *
 */
public class SubgraphAlgorithmMetrics implements SubgraphAssemblyAlgorithmTracker {
	private long timeKmerConstructionStart;
	private int minAnchor;
	private int maxAnchor;
	private long timeSplitOutReferencePathsComplete;
	private int pathsSplit;
	private long timeCalcNonReferenceSubgraphsComplete;
	private int subgraphCount;
	private int startingNodesCounts;
	private int nonreferencePaths;
	private long timeAssemblyNonReferenceContigsComplete;
	private int contigsAssembled;
	private int nodesTraversedDuringContigAssembly;
	private List<BreakendSummary> assemblyBreakends = Lists.newArrayList();
	private long timeShrinkComplete;
	private int shrinkNodesCollapsed;
	private long timeCollapseComplete;
	private int collapseIterations;
	private int pathsCollapsed;
	private int leavesCollapsed;
	private int collapseNodesTraversed;
	private long timeAssemblyStart;
	private long timeAssemblyComplete;
	private long timeGeneratePathGraphComplete;
	private int pathGraphInitialSize;
	private int pathGraphInitialEdges;
	private ProcessingContext processContext;
	private BreakendDirection direction;
	private int referenceIndex;
	private long timeToAssemblyEvidenceComplete;
	private int collapseNodeCountReducedBy;
	private int pathGraphKmers;
	public SubgraphAlgorithmMetrics(ProcessingContext processContext, int referenceIndex, BreakendDirection direction) {
		this(processContext, referenceIndex, direction, System.nanoTime());
	}
	public SubgraphAlgorithmMetrics(ProcessingContext processContext, int referenceIndex, BreakendDirection direction, long constructionStartTime) {
		this.processContext = processContext;
		this.referenceIndex = referenceIndex;
		this.direction = direction;
		this.timeKmerConstructionStart = constructionStartTime;
	}
	/**
	 * Converts the metrics gather into a BED entry
	 * The BED file should contain "gffTags=on" in the track line header 
	 * @return
	 */
	public String toBed() {
		StringBuilder sb = new StringBuilder();
		sb.append(processContext.getDictionary().getSequence(referenceIndex).getSequenceName());
		sb.append('\t');
		sb.append(minAnchor);
		sb.append('\t');
		sb.append(maxAnchor);
		sb.append('\t');
		sb.append(String.format("Name %dms; ", msElapsed(timeKmerConstructionStart, timeAssemblyComplete)));
		sb.append(getGffAttributes());
		sb.append('\t');
		sb.append(msElapsed(timeKmerConstructionStart, timeAssemblyComplete)); 
		sb.append('\t');
		sb.append(direction == BreakendDirection.Forward ? '+' : '-'); // strand
		sb.append('\t');
		sb.append(minAnchor); // thickStart 
		sb.append('\t');
		sb.append(minAnchor); // thickEnd  
		sb.append('\t');
		sb.append(0); // itemRgb  
		sb.append('\t');
		sb.append(assemblyBreakends.size()); // blockCount
		sb.append('\t');
		// blockSizes
		for (BreakendSummary bs : assemblyBreakends) {
			sb.append(bs.end - bs.start + 1);
			sb.append(',');
		}
		sb.append('\t');
		// blockStarts
		for (BreakendSummary bs : assemblyBreakends) {
			sb.append(bs.start - 1); // BED used 0-based indexing
			sb.append(',');
		}
		sb.append('\n');
		return sb.toString();
	}
	/**
	 * Converts to a GFF version 2 
	 * @return
	 */
	public String toGFF2() {
		throw new RuntimeException("NYI");
	}
	private String getGffAttributes() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("Kmers %d; ", pathGraphKmers));
		sb.append(String.format("PathNodes \"%d (%d %d %d)\"; ",
				pathGraphInitialSize - shrinkNodesCollapsed - collapseNodeCountReducedBy + pathsSplit,
				pathGraphInitialSize - shrinkNodesCollapsed,
				-collapseNodeCountReducedBy,
				pathsSplit
				));
		sb.append(String.format("Collapses \"%d (%d steps overs %d iteration)\"; ", pathsCollapsed - leavesCollapsed, collapseNodesTraversed, collapseIterations));
		sb.append(String.format("Subgraphs %d; ", subgraphCount));
		sb.append(String.format("StartingNodes %d; ", startingNodesCounts));
		sb.append(String.format("NonReferencePaths %d; ", nonreferencePaths));
		sb.append(String.format("AssemblyNodeTraversals %d; ", nodesTraversedDuringContigAssembly));
		sb.append(String.format("Times %s; ", toCommaList(toTimeSequenceMs(
				timeKmerConstructionStart,
				timeAssemblyStart,
				timeGeneratePathGraphComplete,
				timeShrinkComplete,
				timeCollapseComplete,
				timeSplitOutReferencePathsComplete,
				timeCalcNonReferenceSubgraphsComplete,
				timeAssemblyNonReferenceContigsComplete,
				timeToAssemblyEvidenceComplete,
				timeAssemblyComplete))));
		return sb.toString();
	}
	private String toCommaList(int... values) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < values.length; i++) {
			sb.append(values[i]);
			if (i < values.length - 1) {
				sb.append(',');
			}
		}
		return sb.toString();
	}
	private int[] toTimeSequenceMs(long... timestamps) {
		int[] out = new int[timestamps.length - 1];
		// strip out zero times
		for (int i = 1; i < timestamps.length; i++) {
			if (timestamps[i] == 0) {
				timestamps[i] = timestamps[i - 1];
			}
			out[i - 1] = msElapsed(timestamps[i - 1], timestamps[i]);
		}
		return out;
	}
	private int msElapsed(long start, long stop) {
		return (int)((stop - start) / 1000000);
	}
	@Override
	public void finalAnchors(int minAnchor, int maxAnchor) {
		this.minAnchor = minAnchor;
		this.maxAnchor = maxAnchor;
	}
	@Override
	public void splitOutReferencePaths(int pathsSplits) {
		this.timeSplitOutReferencePathsComplete = System.nanoTime();
		this.pathsSplit = pathsSplits;
	}
	@Override
	public void calcNonReferenceSubgraphs(
			List<Set<SubgraphPathNode>> subgraphs,
			List<Set<SubgraphPathNode>> startingNodes) {
		this.timeCalcNonReferenceSubgraphsComplete = System.nanoTime();
		this.subgraphCount = subgraphs.size();
		this.startingNodesCounts = 0;
		this.nonreferencePaths = 0;
		for (Set<SubgraphPathNode> s : subgraphs) {
			this.nonreferencePaths += s.size();
			this.startingNodesCounts += startingNodes.size();
		}
	}
	@Override
	public void assemblyNonReferenceContigs(
			List<List<SubgraphPathNode>> assembledPaths,
			List<LinkedList<Long>> assembledKmers, int nodesTraversed) {
		this.timeAssemblyNonReferenceContigsComplete = System.nanoTime();
		this.contigsAssembled = assembledPaths.size();
		this.nodesTraversedDuringContigAssembly = nodesTraversed;
	}
	@Override
	public void toAssemblyEvidence(VariantContextDirectedEvidence assembly) {
		this.timeToAssemblyEvidenceComplete = System.nanoTime();
		this.assemblyBreakends.add(assembly.getBreakendSummary());
	}
	@Override
	public void shrink(int edgesRemoved, int nodesCollapsed) {
		this.timeShrinkComplete = System.nanoTime();
		this.shrinkNodesCollapsed = nodesCollapsed;
	}
	@Override
	public void collapse(int collapseIterations, int pathsCollapsed,
			int leavesCollapsed, int nodesTraversed, int nodeCountReducedBy) {
		this.timeCollapseComplete = System.nanoTime();
		this.collapseIterations = collapseIterations;
		this.pathsCollapsed = pathsCollapsed;
		this.leavesCollapsed = leavesCollapsed;
		this.collapseNodesTraversed = nodesTraversed;
		this.collapseNodeCountReducedBy = nodeCountReducedBy;
	}
	@Override
	public void assemblyStarted() {
		this.timeAssemblyStart = System.nanoTime();
	}
	@Override
	public void assemblyComplete() {
		this.timeAssemblyComplete = System.nanoTime();
	}
	@Override
	public void generatePathGraph(int kmers, int nodes, int edges) {
		this.timeGeneratePathGraphComplete = System.nanoTime();
		this.pathGraphKmers = kmers;
		this.pathGraphInitialSize = nodes;
		this.pathGraphInitialEdges = edges;
	}
	@Override
	public int getReferenceIndex() {
		return referenceIndex;
	}
	@Override
	public long getStartAnchorPosition() {
		return minAnchor;
	}
	@Override
	public long getEndAnchorPosition() {
		return maxAnchor;
	}
}
