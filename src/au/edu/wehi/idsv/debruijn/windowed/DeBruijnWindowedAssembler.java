package au.edu.wehi.idsv.debruijn.windowed;

import htsjdk.samtools.SAMRecord;

import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceEndCoordinateComparator;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
import au.edu.wehi.idsv.ReadEvidenceAssemblerUtil;
import au.edu.wehi.idsv.SAMRecordUtil;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.sam.AnomolousReadAssembly;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Generates local breakpoint read de bruijn graph assemblies of SV-supporting reads
 * in a single coordinate-sorted pass over the read evidence.
 * 
 * An assembly is generated for each non-reference contig.
 * 
 * @author Daniel Cameron
 *
 */
public class DeBruijnWindowedAssembler implements ReadEvidenceAssembler {
	private final ProcessingContext processContext;
	private final int k;
	/**
	 * Array size is 2: forward then backward graphs
	 */
	private DeBruijnReadGraph[] graphs;
	private int currentReferenceIndex = -1;
	public DeBruijnWindowedAssembler(
			ProcessingContext processContext,
			int kmerSize) {
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
		for (DeBruijnReadGraph g : graphs) {
			if (g.getDirection() == evidence.getBreakendSummary().direction) {
				g.addEvidence(evidence);
			}
		}
		int startpos = evidence.getBreakendSummary().start;
		// make sure we have enough margin that we don't assemble a variant before all the
		// evidence for it is available
		int shouldBeCompletedPos = startpos - 3 * processContext.getMetrics().getMaxFragmentSize();
		it = Iterables.concat(Iterables.mergeSorted(
				graphs[0].assembleContigsBefore(shouldBeCompletedPos),
				graphs[1].assembleContigsBefore(shouldBeCompletedPos)
				));
		
		// TODO: sort by position
		return it;
	}
	@Override
	public Iterable<VariantContextDirectedBreakpoint> endOfEvidence() {
		return assembleAll();
	}
	private void init(int referenceIndex) {
		currentReferenceIndex = referenceIndex;
		graphs = new DeBruijnReadGraph[] {
			new DeBruijnReadGraph(processContext, referenceIndex, k, BreakendDirection.Forward),
			new DeBruijnReadGraph(processContext, referenceIndex, k, BreakendDirection.Backward),
		};
	}
	private Iterable<VariantContextDirectedBreakpoint> assembleAll() {
		if (graphs == null) return ImmutableList.of();
		return Iterables.mergeSorted(graphs[0].assembleContigsBefore(Integer.MAX_VALUE), graphs[1].assembleContigsBefore(Integer.MAX_VALUE));
		// done with these graphs
		graphs = null;
	}
}
