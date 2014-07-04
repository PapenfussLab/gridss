package au.edu.wehi.idsv.debruijn.subgraph;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

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
	private final int k;
	/**
	 * Array size is 2: forward then backward graphs
	 */
	private DeBruijnReadGraph fgraph;
	private DeBruijnReadGraph bgraph;
	private int currentReferenceIndex = -1;
	public DeBruijnSubgraphAssembler(ProcessingContext processContext, int kmerSize) {
		this.processContext = processContext;
		this.k = kmerSize;
	}
	@Override
	public Iterable<VariantContextDirectedBreakpoint> addEvidence(DirectedEvidence evidence) {
		Iterable<VariantContextDirectedBreakpoint> it = ImmutableList.of();
		if (evidence.getBreakendSummary().referenceIndex != currentReferenceIndex) {
			it = assembleAll();
			init(evidence.getBreakendSummary().referenceIndex);
		}
		if (evidence.getBreakendSummary().direction == BreakendDirection.Forward) {
			fgraph.addEvidence(evidence);
		} else {
			bgraph.addEvidence(evidence);
		}
		int startpos = evidence.getBreakendSummary().start;
		// make sure we have enough margin that we don't assemble a variant before all the
		// evidence for it is available
		int shouldBeCompletedPos = startpos - 3 * processContext.getMetrics().getMaxFragmentSize();
		it = Iterables.concat(it, assembleBefore(shouldBeCompletedPos));
		return it;
	}
	@Override
	public Iterable<VariantContextDirectedBreakpoint> endOfEvidence() {
		return assembleAll();
	}
	private void init(int referenceIndex) {
		currentReferenceIndex = referenceIndex;
		fgraph = new DeBruijnReadGraph(processContext, referenceIndex, k, BreakendDirection.Forward);
		bgraph = new DeBruijnReadGraph(processContext, referenceIndex, k, BreakendDirection.Backward);
	}
	private Iterable<VariantContextDirectedBreakpoint> assembleAll() {
		return assembleBefore(Integer.MAX_VALUE);
	}
	private Iterable<VariantContextDirectedBreakpoint> assembleBefore(int position) {
		if (currentReferenceIndex < 0) return ImmutableList.of();
		Iterable<VariantContextDirectedBreakpoint> it = Iterables.mergeSorted(ImmutableList.of(
				fgraph.assembleContigsBefore(position),
				bgraph.assembleContigsBefore(position)),
				VariantContextDirectedBreakpoint.ByLocation);
		fgraph.removeBefore(position);
		bgraph.removeBefore(position);
		return it;
	}
}
