package au.edu.wehi.idsv.debruijn.subgraph;

import java.io.File;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.visualisation.DeBruijnSubgraphGexfExporter;

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
	private final AssemblyParameters parameters;
	private DeBruijnReadGraph fgraph;
	private DeBruijnReadGraph bgraph;
	private int currentReferenceIndex = -1;
	private int processStep = 0;
	public DeBruijnSubgraphAssembler(ProcessingContext processContext, AssemblyEvidenceSource source) {
		this.processContext = processContext;
		this.source = source;
		this.parameters = processContext.getAssemblyParameters();
	}
	@Override
	public Iterable<VariantContextDirectedEvidence> addEvidence(DirectedEvidence evidence) {
		Iterable<VariantContextDirectedEvidence> it = ImmutableList.of();
		if (evidence.getBreakendSummary().referenceIndex != currentReferenceIndex) {
			it = assembleAll();
			init(evidence.getBreakendSummary().referenceIndex);
		}
		// Assemble old evidence that couldn't have anything to do with us
		int startpos = evidence.getBreakendSummary().start;
		int shouldBeCompletedPos = (int)(startpos - processContext.getAssemblyParameters().subgraphAssemblyMargin * source.getMetrics().getMaxFragmentSize());
		it = Iterables.concat(it, assembleBefore(shouldBeCompletedPos));
		if (DeBruijnGraphBase.PERFORM_EXPENSIVE_SANITY_CHECKS) {
			assert(fgraph.sanityCheckSubgraphs(shouldBeCompletedPos, startpos + source.getMetrics().getMaxFragmentSize()));
			assert(bgraph.sanityCheckSubgraphs(shouldBeCompletedPos, startpos + source.getMetrics().getMaxFragmentSize()));
		}
		startingNextProcessingStep();
		if (evidence.getBreakendSummary().direction == BreakendDirection.Forward) {
			fgraph.addEvidence(evidence);
		} else {
			bgraph.addEvidence(evidence);
		}
		return it;
	}
	@Override
	public Iterable<VariantContextDirectedEvidence> endOfEvidence() {
		return assembleAll();
	}
	private void init(int referenceIndex) {
		processStep = 0;
		currentReferenceIndex = referenceIndex;
		fgraph = new DeBruijnReadGraph(processContext, source, referenceIndex, BreakendDirection.Forward, parameters);
		bgraph = new DeBruijnReadGraph(processContext, source, referenceIndex, BreakendDirection.Backward, parameters);
		if (parameters.debruijnGraphVisualisationDirectory != null) {
			fgraph.setGraphExporter(new DeBruijnSubgraphGexfExporter(parameters.k));
			bgraph.setGraphExporter(new DeBruijnSubgraphGexfExporter(parameters.k));
		}
	}
	private Iterable<VariantContextDirectedEvidence> assembleAll() {
		Iterable<VariantContextDirectedEvidence> assemblies = assembleBefore(Integer.MAX_VALUE);
		File exportDir = parameters.debruijnGraphVisualisationDirectory;
		if (exportDir != null) {
			exportDir.mkdir();
			if (fgraph != null) {
				fgraph.getGraphExporter().saveTo(new File(exportDir, String.format("debruijn.kmers.forward.%s.gexf", processContext.getDictionary().getSequence(currentReferenceIndex).getSequenceName())));
			}
			if (bgraph != null) {
				bgraph.getGraphExporter().saveTo(new File(exportDir, String.format("debruijn.kmers.backward.%s.gexf", processContext.getDictionary().getSequence(currentReferenceIndex).getSequenceName())));
			}
		}
		fgraph = null;
		bgraph = null;
		return assemblies;
	}
	private Iterable<VariantContextDirectedEvidence> assembleBefore(int position) {
		startingNextProcessingStep();
		if (currentReferenceIndex < 0) return ImmutableList.of();
		Iterable<VariantContextDirectedEvidence> it = Iterables.mergeSorted(ImmutableList.of(
				fgraph.assembleContigsBefore(position),
				bgraph.assembleContigsBefore(position)),
				VariantContextDirectedEvidence.ByLocationStart);
		fgraph.removeBefore(position);
		bgraph.removeBefore(position);
		return it;
	}
	/**
	 * Indicates that the next processing step has started
	 */
	private void startingNextProcessingStep() {
		processStep++;
		if (fgraph != null && fgraph.getGraphExporter() != null) fgraph.getGraphExporter().setTime(processStep);
		if (bgraph != null && bgraph.getGraphExporter() != null) bgraph.getGraphExporter().setTime(processStep);
	}
}
