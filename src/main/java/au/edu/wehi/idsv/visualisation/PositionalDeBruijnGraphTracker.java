package au.edu.wehi.idsv.visualisation;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.positional.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.*;

/**
 * Tracks information associated with positional de Bruijn graph calling
 * @author Daniel Cameron
 *
 */
public class PositionalDeBruijnGraphTracker implements Closeable {
	private Log log = Log.getInstance(PositionalDeBruijnGraphTracker.class);
	public static class ContigStats {
		public int contigStartPosition;
		public int contigNodes;
		public int startAnchorNodes;
		public int endAnchorNodes;
		public int truncatedNodes;
		public static String header() {
			return "contigStartPosition,contigNodeSize,contigStartAnchorNodeSize,contigEndAnchorNodeSize,contigTruncatedNodeSize";
		}
		@Override
		public String toString() {
			return String.format("%d,%d,%d,%d,%d", contigStartPosition,contigNodes,startAnchorNodes,endAnchorNodes,truncatedNodes);
		}
	}
	public static class MemoizationStats {
		public int nodes;
		public int removed;
		public int pathsRemoved;
		public int descendentPathsRemoved;
		public int pathsReset;
		public int pathsRestarted;
		public static String header() {
			return "memoizedSize,memoizedRemovalSize,memoizedPathsRemovalSize,descendentPathsRemovalSize,memoizedPathsReactivateSize,memoizedPathsRestartSize";
		}
		@Override
		public String toString() {
			return String.format("%d,%d,%d,%d,%d,%d", nodes, removed, pathsRemoved, descendentPathsRemoved, pathsReset, pathsRestarted);
		}
	}
	private BufferedWriter writer;
	private SupportNodeIterator support;
	private AggregateNodeIterator aggregate;
	private PathNodeIterator pathNode;
	private CollapseIterator collapse;
	private PathSimplificationIterator simplify;
	private EvidenceTracker tracker;
	private NonReferenceContigAssembler assembler;
	private long lastTime = System.nanoTime();
	public PositionalDeBruijnGraphTracker(
			File file,
			SupportNodeIterator support,
			AggregateNodeIterator aggregate,
			PathNodeIterator pathNode,
			CollapseIterator collapse,
			PathSimplificationIterator simplify,
			EvidenceTracker tracker,
			NonReferenceContigAssembler assembler) throws IOException {
		this.writer = new BufferedWriter(new FileWriter(file));
		this.support = support;
		this.pathNode = pathNode;
		this.aggregate = aggregate;
		this.collapse = collapse;
		this.simplify = simplify;
		this.tracker = tracker;
		this.assembler = assembler;
	}
	public void writeHeader() throws IOException {
		if (writer == null) return;
		writer.write("nsElapsedTime");
		writer.write(",supportPosition,aggregatePosition,pathNodePosition,collapsePosition,simplifyPosition,assemblerPosition,assemblerFirstPosition");
		writer.write(",supportConsumed,aggregateConsumed,pathNodeConsumed,collapseConsumed,simplifyConsumed,assemblerConsumed,trackerConsumed");
		writer.write(",trackerActive");
		writer.write(",supportProcessedSize");
		writer.write(",aggregateProcessedSize,aggregateQueueSize,aggregateActiveSize");
		writer.write(",pathNodeProcessedSize,pathNodeActiveSize,pathNodeEdgeLookupSize,pathNodePathLookupSize");
		writer.write(",collapseProcessedSize,collapseUnprocessedSize,collapseTraversalCount,collapsedBranchCount,collapsedLeafCount");
		writer.write(",simplifyProcessedSize,simplifyLookupSize,simplifyUnprocessedSize,simplifiedCount");
		writer.write(",trackerLookupSize");
		writer.write(",contigFrontierSize,contigMemoizedSize");
		writer.write(",assemblyActiveSize");
		writer.write(",");
		writer.write(ContigStats.header());
		writer.write(",");
		writer.write(MemoizationStats.header());
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			writer.write(",aggregateKmerMaxActive,aggregateActiveNodes,pathNodeEdgeMaxActive,pathNodePathMaxActive,trackerMaxKmerSupport,assemblyMaxActive,trackerLookupSize");
		}
		writer.write('\n');
	}
	public void trackAssembly(MemoizedContigCaller caller) {
		if (writer == null) return;
		long currentTime = System.nanoTime();
		long deltaTime = currentTime - lastTime;
		try {
			writer.write(Long.toString(deltaTime));
			writer.write(',');
			writer.write(Integer.toString(support.tracking_inputPosition()));
			writer.write(',');
			writer.write(Integer.toString(aggregate.tracking_inputPosition()));
			writer.write(',');
			writer.write(Integer.toString(pathNode.tracking_inputPosition()));
			writer.write(',');
			if (collapse != null) writer.write(Integer.toString(collapse.tracking_inputPosition()));
			writer.write(',');
			if (simplify != null) writer.write(Integer.toString(simplify.tracking_inputPosition()));
			writer.write(',');
			writer.write(Integer.toString(assembler.tracking_firstPosition()));
			writer.write(',');
			writer.write(Integer.toString(assembler.tracking_inputPosition()));
			writer.write(',');
			writer.write(Long.toString(support.tracking_underlyingConsumed()));
			writer.write(',');
			writer.write(Long.toString(aggregate.tracking_underlyingConsumed()));
			writer.write(',');
			writer.write(Long.toString(pathNode.tracking_underlyingConsumed()));
			writer.write(',');
			if (collapse != null) writer.write(Long.toString(collapse.tracking_underlyingConsumed()));
			writer.write(',');
			if (simplify != null) writer.write(Long.toString(simplify.tracking_underlyingConsumed()));
			writer.write(',');
			writer.write(Long.toString(assembler.tracking_underlyingConsumed()));
			writer.write(',');
			writer.write(Long.toString(tracker.tracking_evidenceTotal()));
			writer.write(',');
			writer.write(Long.toString(tracker.tracking_evidenceActive()));
			writer.write(',');
			writer.write(Integer.toString(support.tracking_processedSize()));
			writer.write(',');
			writer.write(Integer.toString(aggregate.tracking_processedSize()));
			writer.write(',');
			writer.write(Integer.toString(aggregate.tracking_aggregatorQueueSize()));
			writer.write(',');
			writer.write(Integer.toString(aggregate.tracking_kmerCount()));
			writer.write(',');
			writer.write(Integer.toString(pathNode.tracking_processedSize()));
			writer.write(',');
			writer.write(Integer.toString(pathNode.tracking_activeSize()));
			writer.write(',');
			writer.write(Integer.toString(pathNode.tracking_edgeLookupSize()));
			writer.write(',');
			writer.write(Integer.toString(pathNode.tracking_pathNodeEdgeLookupSize()));
			writer.write(',');
			if (collapse != null) writer.write(Integer.toString(collapse.tracking_processedSize()));
			writer.write(',');
			if (collapse != null) writer.write(Integer.toString(collapse.tracking_unprocessedSize()));
			writer.write(',');
			if (collapse != null) writer.write(Long.toString(collapse.tracking_traversalCount()));
			writer.write(',');
			if (collapse != null) writer.write(Long.toString(collapse.tracking_branchCollapseCount()));
			writer.write(',');
			if (collapse != null) writer.write(Long.toString(collapse.tracking_leafCollapseCount()));
			writer.write(',');
			if (simplify != null) writer.write(Integer.toString(simplify.tracking_processedSize()));
			writer.write(',');
			if (simplify != null) writer.write(Integer.toString(simplify.tracking_lookupSize()));
			writer.write(',');
			if (simplify != null) writer.write(Integer.toString(simplify.tracking_unprocessedSize()));
			writer.write(',');
			if (simplify != null) writer.write(Long.toString(simplify.tracking_simplifiedCount()));
			writer.write(',');
			writer.write(Integer.toString(tracker.tracking_kmerCount()));
			writer.write(',');
			writer.write(Integer.toString(caller.tracking_frontierSize()));
			writer.write(',');
			writer.write(Integer.toString(caller.tracking_memoizedNodeCount()));
			writer.write(',');
			writer.write(Integer.toString(assembler.tracking_activeNodes()));
			writer.write(',');
			writer.write(assembler.tracking_lastContig().toString());
			writer.write(',');
			writer.write(caller.tracking_lastRemoval().toString());
			if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
				writer.write(',');
				writer.write(Integer.toString(aggregate.tracking_aggregatorKmerMaxActiveNodeCount()));
				writer.write(',');
				writer.write(Integer.toString(aggregate.tracking_aggregatorActiveNodeCount()));
				writer.write(',');
				writer.write(Integer.toString(pathNode.tracking_pathNodeEdgeLookupMaxKmerNodeCount()));
				writer.write(',');
				writer.write(Integer.toString(pathNode.tracking_edgeLookupMaxKmerNodeCount()));
				writer.write(',');
				writer.write(Integer.toString(tracker.tracking_maxKmerSupportNodesCount()));
				writer.write(',');
				writer.write(Integer.toString(assembler.tracking_maxKmerActiveNodeCount()));
				writer.write(',');
				writer.write(Integer.toString(tracker.tracking_supportNodeCount()));	
			}
			writer.write('\n');
		} catch (IOException e) {
			if (log != null) log.error(e);
			log = null;
		} finally {
			lastTime = currentTime;
		}
	}
	@Override
	public void close() throws IOException {
		if (writer != null) writer.flush();
		CloserUtil.close(writer);
		writer = null;
		support = null;
		aggregate = null;
		tracker = null;
		assembler = null;
	}
}
