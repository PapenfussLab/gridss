package au.edu.wehi.idsv.debruijn.anchored;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
import au.edu.wehi.idsv.RealignedRemoteSoftClipEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;

/**
 * Generates local breakpoint read de bruijn graph assemblies of SV-supporting reads
 * in a single coordinate-sorted pass over the read evidence.
 * 
 * An anchored assembly for each breakend direction is assembled for each position.
 * 
 * @author Daniel Cameron
 *
 */
public class DeBruijnAnchoredAssembler implements ReadEvidenceAssembler {
	private final DeBruijnAnchoredGraph graphf;
	private final DeBruijnAnchoredGraph graphb;
	private final PriorityQueue<DirectedEvidence> activef = new PriorityQueue<DirectedEvidence>(1024, DirectedEvidenceOrder.ByEndStart);
	private final PriorityQueue<DirectedEvidence> activeb = new PriorityQueue<DirectedEvidence>(1024, DirectedEvidenceOrder.ByEndStart);
	private int currentReferenceIndex = -1;
	private int currentPosition = -1;
	public DeBruijnAnchoredAssembler(ProcessingContext processContext, AssemblyEvidenceSource source) {
		this.graphf = new DeBruijnAnchoredGraph(processContext, source, processContext.getAssemblyParameters().k, BreakendDirection.Forward);
		this.graphb = new DeBruijnAnchoredGraph(processContext, source, processContext.getAssemblyParameters().k, BreakendDirection.Backward);
	}
	private AssemblyEvidence createAssemblyBreakend(DeBruijnAnchoredGraph graph, BreakendDirection direction) {
		AssemblyEvidence variant = graph.assemble(currentReferenceIndex, currentPosition);
		return variant;
	}
	@Override
	public Iterable<AssemblyEvidence> addEvidence(DirectedEvidence evidence) {
		if (evidence == null) return Collections.emptyList();
		BreakendSummary location = evidence.getBreakendSummary();
		List<AssemblyEvidence> result = processUpToExcluding(location.referenceIndex, location.start);
		// add evidence
		if (location.direction == BreakendDirection.Forward) {
			addToGraph(graphf, activef, evidence);
		} else {
			addToGraph(graphb, activeb, evidence);
		}
		return result;
	}
	/**
	 * Adds the given evidence to the given graph
	 * @param graph graph to add evidence to
	 * @param active reads current active in the given graph
	 * @param evidence read evidence to add
	 */
	private static void addToGraph(DeBruijnAnchoredGraph graph, PriorityQueue<DirectedEvidence> active, DirectedEvidence evidence) {
		if (evidence instanceof RealignedRemoteSoftClipEvidence) {
			// ignore remote soft clips - they get assembled on the other side
		} else if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair nrrp = (NonReferenceReadPair)evidence;
			graph.addEvidence(nrrp);
			active.add(evidence);
		} else if (evidence instanceof SoftClipEvidence) {
			SoftClipEvidence sc = (SoftClipEvidence)evidence;
			graph.addEvidence(sc);
			active.add(evidence);
		}
	}
	private static void flushUpToExcluding(DeBruijnAnchoredGraph graph, PriorityQueue<DirectedEvidence> active, int referenceIndex, int position) {
		while (!active.isEmpty() &&
				(active.peek().getBreakendSummary().referenceIndex < referenceIndex || 
				(active.peek().getBreakendSummary().referenceIndex == referenceIndex && active.peek().getBreakendSummary().end < position))) {
			DirectedEvidence evidence = active.poll();
			if (evidence instanceof NonReferenceReadPair) {
				graph.removeEvidence((NonReferenceReadPair)evidence);
			} else if (evidence instanceof SoftClipEvidence) {
				graph.removeEvidence((SoftClipEvidence)evidence);
			} else {
				throw new IllegalStateException(String.format("Sanity check failure: unhandled evidence of type %s present in de bruijn graph", evidence.getClass()));
			}
		}
	}
	@Override
	public Iterable<AssemblyEvidence> endOfEvidence() {
		return processUpToExcluding(currentReferenceIndex + 1, 0);
	}
	/**
	 * Processes all positions before the given position
	 * @param referenceIndex reference index of position to process until
	 * @param position position to process until
	 * @return add assemblies generated up to the given position
	 */
	private List<AssemblyEvidence> processUpToExcluding(int referenceIndex, int position) {
		List<AssemblyEvidence> result = new ArrayList<>();
		// keep going till we flush both de bruijn graphs
		// or we reach the stop position
		while ((!activef.isEmpty() || !activeb.isEmpty()) &&
				!(currentReferenceIndex == referenceIndex && currentPosition == position)) {
			AssemblyEvidence assembly;
			
			assembly = createAssemblyBreakend(graphf, BreakendDirection.Forward);
			if (assembly != null) {
				result.add(assembly);
			}
			assembly = createAssemblyBreakend(graphb, BreakendDirection.Backward);
			if (assembly != null) {
				result.add(assembly);
			}
			currentPosition++;
			flushUpToExcluding(graphf, activef, currentReferenceIndex, currentPosition);
			flushUpToExcluding(graphb, activeb, currentReferenceIndex, currentPosition);
		}
		currentReferenceIndex = referenceIndex;
		currentPosition = position;
		return result;
	}
	@Override
	public String getStateSummaryMetrics() {
		return "";
	}
}
