package au.edu.wehi.idsv.debruijn.subgraph;

import java.io.File;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
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
	private SubgraphAssemblyAlgorithmTrackerBEDWriter currentTracker = null;
	private int processStep = 0;
	public DeBruijnSubgraphAssembler(ProcessingContext processContext, AssemblyEvidenceSource source) {
		this.processContext = processContext;
		this.source = source;
	}
	@Override
	public Iterable<AssemblyEvidence> addEvidence(DirectedEvidence evidence) {
		if (evidence.getBreakendSummary() == null) throw new IllegalArgumentException("Invalid evidence");
		Iterable<AssemblyEvidence> it = ImmutableList.of();
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
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(graph.sanityCheckSubgraphs());
		}
		return it;
	}
	@Override
	public Iterable<AssemblyEvidence> endOfEvidence() {
		return assembleAll();
	}
	private void init(int referenceIndex) {
		processStep = 0;
		currentReferenceIndex = referenceIndex;
		if (currentTracker != null) {
			currentTracker.close();
		}
		if (processContext.getAssemblyParameters().trackAlgorithmProgress) {
			currentTracker = new SubgraphAssemblyAlgorithmTrackerBEDWriter(
					(int)(processContext.getAssemblyParameters().subgraphAssemblyMargin * source.getMaxConcordantFragmentSize()),
				new File(processContext.getAssemblyParameters().debruijnGraphVisualisationDirectory,
					String.format("debruijn.assembly.metrics.%s.bed", processContext.getDictionary().getSequence(currentReferenceIndex).getSequenceName())));
		}
		graph = new DeBruijnReadGraph(processContext, source, referenceIndex, processContext.getAssemblyParameters(), currentTracker);
		if (processContext.getAssemblyParameters().debruijnGraphVisualisationDirectory != null && processContext.getAssemblyParameters().visualiseAll) {
			graph.setGraphExporter(new DeBruijnSubgraphGexfExporter(processContext.getAssemblyParameters().k));
		
		}
	}
	private Iterable<AssemblyEvidence> assembleAll() {
		Iterable<AssemblyEvidence> assemblies = assembleBefore(Long.MAX_VALUE);
		File exportDir = processContext.getAssemblyParameters().debruijnGraphVisualisationDirectory;
		if (exportDir != null && processContext.getAssemblyParameters().visualiseAll) {
			exportDir.mkdir();
			if (graph != null) {
				graph.getGraphExporter().saveTo(new File(exportDir, String.format("debruijn.kmers.%s.gexf", processContext.getDictionary().getSequence(currentReferenceIndex).getSequenceName())));
			}
		}
		graph = null;
		if (currentTracker != null) {
			currentTracker.close();
			currentTracker = null;
		}
		return assemblies;
	}
	private Iterable<AssemblyEvidence> assembleBefore(long position) {
		startingNextProcessingStep();
		if (currentReferenceIndex < 0) return ImmutableList.of();
		Iterable<AssemblyEvidence> it = graph.assembleContigsBefore(position);
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
