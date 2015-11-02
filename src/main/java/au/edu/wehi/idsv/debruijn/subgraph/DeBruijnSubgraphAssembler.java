package au.edu.wehi.idsv.debruijn.subgraph;

import java.io.File;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNode;
import au.edu.wehi.idsv.visualisation.DeBruijnSubgraphGexfExporter;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTrackerBEDWriter;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;

/**
 * Generates local breakpoint read de bruijn graph assemblies of SV-supporting reads
 * in a single coordinate-sorted pass over the read evidence.
 * 
 * Assembly is performed for each subgraph, with each non-reference kmer contributing
 * to at most one assembly
 * 
 * @author Daniel Cameron
 *
 */
public class DeBruijnSubgraphAssembler implements ReadEvidenceAssembler {
	private final ProcessingContext processContext;
	private final AssemblyEvidenceSource source;
	private DeBruijnReadGraph graph;
	private int currentReferenceIndex = -1;
	private SubgraphAssemblyAlgorithmTrackerBEDWriter<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> currentTracker = null;
	private int processStep = 0;
	public DeBruijnSubgraphAssembler(ProcessingContext processContext, AssemblyEvidenceSource source) {
		this.processContext = processContext;
		this.source = source;
	}
	@Override
	public Iterable<SAMRecordAssemblyEvidence> addEvidence(DirectedEvidence evidence) {
		if (evidence.getBreakendSummary() == null) throw new IllegalArgumentException("Invalid evidence");
		Iterable<SAMRecordAssemblyEvidence> it = ImmutableList.of();
		if (evidence.getBreakendSummary().referenceIndex != currentReferenceIndex) {
			it = assembleAll();
			init(evidence.getBreakendSummary().referenceIndex);
		}
		// Assemble old evidence that couldn't have anything to do with us
		assert(evidence.getBreakendSummary().referenceIndex == currentReferenceIndex);
		long startpos = processContext.getLinear().getLinearCoordinate(currentReferenceIndex, evidence.getBreakendSummary().start);
		long shouldBeCompletedPos = startpos - source.getAssemblyEvidenceWindowSize();
		it = Iterables.concat(it, assembleBefore(shouldBeCompletedPos));
		startingNextProcessingStep();
		graph.addEvidence(evidence);
		if (Defaults.SANITY_CHECK_DE_BRUIJN) {
			assert(graph.sanityCheckSubgraphs());
		}
		return it;
	}
	@Override
	public Iterable<SAMRecordAssemblyEvidence> endOfEvidence() {
		return assembleAll();
	}
	private void init(int referenceIndex) {
		processStep = 0;
		currentReferenceIndex = referenceIndex;
		if (currentTracker != null) {
			currentTracker.close();
		}
		if (processContext.getConfig().getVisualisation().assemblyProgress) {
			currentTracker = new SubgraphAssemblyAlgorithmTrackerBEDWriter<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>>(
					(int)(processContext.getAssemblyParameters().subgraph.assemblyMargin * source.getMaxConcordantFragmentSize()),
				new File(processContext.getConfig().getVisualisation().directory,
					String.format("debruijn.assembly.metrics.%s.bed", processContext.getDictionary().getSequence(currentReferenceIndex).getSequenceName())));
		}
		graph = new DeBruijnReadGraph(source, referenceIndex, currentTracker);
		if (processContext.getConfig().getVisualisation().assembly) {
			graph.setGraphExporter(new DeBruijnSubgraphGexfExporter(processContext.getAssemblyParameters().k));
		
		}
	}
	private Iterable<SAMRecordAssemblyEvidence> assembleAll() {
		Iterable<SAMRecordAssemblyEvidence> assemblies = assembleBefore(Long.MAX_VALUE);
		if (processContext.getConfig().getVisualisation().assembly) {
			if (graph != null) {
				graph.getGraphExporter().saveTo(new File(processContext.getConfig().getVisualisation().directory, String.format("debruijn.kmers.%s.gexf", processContext.getDictionary().getSequence(currentReferenceIndex).getSequenceName())));
			}
		}
		graph = null;
		if (currentTracker != null) {
			currentTracker.close();
			currentTracker = null;
		}
		return assemblies;
	}
	private Iterable<SAMRecordAssemblyEvidence> assembleBefore(long position) {
		startingNextProcessingStep();
		if (currentReferenceIndex < 0) return ImmutableList.of();
		Iterable<SAMRecordAssemblyEvidence> it = graph.assembleContigsBefore(position);
		graph.removeBefore(position);
		return it;
	}
	/**
	 * Indicates that the next processing step has started
	 */
	private void startingNextProcessingStep() {
		processStep++;
		if (graph != null && graph.getGraphExporter() != null) graph.getGraphExporter().setTime(processStep);
	}
	@Override
	public String getStateSummaryMetrics() {
		return graph.getStateSummaryMetrics();
	}
}
