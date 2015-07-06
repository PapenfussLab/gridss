package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.debruijn.positional.AggregateNodeIterator;
import au.edu.wehi.idsv.debruijn.positional.BestNonReferenceContigCaller;
import au.edu.wehi.idsv.debruijn.positional.EvidenceTracker;
import au.edu.wehi.idsv.debruijn.positional.KmerPathSubnode;
import au.edu.wehi.idsv.debruijn.positional.NonReferenceContigAssembler;
import au.edu.wehi.idsv.debruijn.positional.PathCollapseIterator;
import au.edu.wehi.idsv.debruijn.positional.PathNodeIterator;
import au.edu.wehi.idsv.debruijn.positional.PathSimplificationIterator;
import au.edu.wehi.idsv.debruijn.positional.SupportNodeIterator;

/**
 * Tracks information associated with positional de Bruijn graph calling
 * @author cameron.d
 *
 */
public class PositionalDeBruijnGraphTracker implements Closeable {
	private Log log = Log.getInstance(PositionalDeBruijnGraphTracker.class);
	private BufferedWriter writer;
	private SupportNodeIterator support;
	private AggregateNodeIterator aggregate;
	private PathNodeIterator pathNode;
	private PathCollapseIterator collapse;
	private PathSimplificationIterator simplify;
	private EvidenceTracker tracker;
	private NonReferenceContigAssembler assembler;
	public PositionalDeBruijnGraphTracker(
			File file,
			SupportNodeIterator support,
			AggregateNodeIterator aggregate,
			PathNodeIterator pathNode,
			PathCollapseIterator collapse,
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
		writer.write("supportPosition,aggregatePosition,pathNodePosition,collapsePosition,simplifyPosition,firstContigPosition,calledContigPosition,assemblerPosition,assemblerFirstPosition");
		writer.write(",supportConsumed,aggregateConsumed,pathNodeConsumed,collapseConsumed,simplifyConsumed,contigConsumed,assemblerConsumed,trackerConsumed");
		writer.write(",supportProcessedSize");
		writer.write(",aggregateProcessedSize,aggregateQueueSize,aggregateActiveSize");
		writer.write(",pathNodeProcessedSize,pathNodeActiveSize,pathNodeEdgeLookupSize,pathNodePathLookupSize");
		writer.write(",collapseProcessedSize,collapseUnprocessedSize,collapseTraversalCount,collapsedBranchCount,collapsedLeafCount");
		writer.write(",simplifyProcessedSize,simplifyLookupSize,simplifyUnprocessedSize,simplifiedCount");
		writer.write(",trackerLookupSize");
		writer.write(",contigSize,contigFrontierSize,contigMemoizedSize,contigUnprocessedSize");
		writer.write(",assemblyActiveSize");
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			writer.write(",aggregateKmerMaxActive,aggregateActiveNodes,pathNodeEdgeMaxActive,pathNodePathMaxActive,trackerMaxKmerSupport,assemblyMaxActive,trackerLookupSize");
		}
		writer.write('\n');
	}
	public void trackAssembly(BestNonReferenceContigCaller caller) {
		if (writer == null) return;
		try {
			ArrayDeque<KmerPathSubnode> called = caller.bestContig();
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
			writer.write(Integer.toString(caller.tracking_contigFirstPosition()));
			writer.write(',');
			if (called != null) writer.write(Integer.toString(called.getFirst().firstStart()));
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
			writer.write(Long.toString(caller.tracking_underlyingConsumed()));
			writer.write(',');
			writer.write(Long.toString(assembler.tracking_underlyingConsumed()));
			writer.write(',');
			writer.write(Long.toString(tracker.tracking_underlyingConsumed()));
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
			writer.write(Integer.toString(caller.tracking_contigCount()));
			writer.write(',');
			writer.write(Integer.toString(caller.tracking_frontierSize()));
			writer.write(',');
			writer.write(Integer.toString(caller.tracking_memoizedNodeCount()));
			writer.write(',');
			writer.write(Integer.toString(caller.tracking_unprocessedStartNodeCount()));
			writer.write(',');
			writer.write(Integer.toString(assembler.tracking_activeNodes()));
			
			if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
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
