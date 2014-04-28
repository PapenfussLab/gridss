package au.edu.wehi.socrates.debruijn;

import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import au.edu.wehi.socrates.BreakpointDirection;
import au.edu.wehi.socrates.BreakpointLocation;
import au.edu.wehi.socrates.DirectedBreakpointAssembly;
import au.edu.wehi.socrates.DirectedEvidence;
import au.edu.wehi.socrates.DirectedEvidenceEndCoordinateComparator;
import au.edu.wehi.socrates.LinearGenomicCoordinate;
import au.edu.wehi.socrates.NonReferenceReadPair;
import au.edu.wehi.socrates.ReadEvidenceAssembler;
import au.edu.wehi.socrates.SoftClipEvidence;
import au.edu.wehi.socrates.sam.AnomolousReadAssembly;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Queues;

/**
 * Generates local breakpoint read de bruijn graph assemblies of SV-supporting reads
 * in a single coordinate-sorted pass over the read evidence
 * @author Daniel Cameron
 *
 */
public class DeBruijnAssembler implements ReadEvidenceAssembler {
	public static final String ASSEMBLER_NAME = "debruijn";
	private final SAMSequenceDictionary dictionary;
	private final ReferenceSequenceFile reference;
	private final DeBruijnReadGraph graphf;
	private final DeBruijnReadGraph graphb;
	private final PriorityQueue<DirectedEvidence> activef = new PriorityQueue<DirectedEvidence>(1024, new DirectedEvidenceEndCoordinateComparator());
	private final PriorityQueue<DirectedEvidence> activeb = new PriorityQueue<DirectedEvidence>(1024, new DirectedEvidenceEndCoordinateComparator());
	private int currentReferenceIndex = -1;
	private int currentPosition = -1;
	public DeBruijnAssembler(
			SAMSequenceDictionary dictionary,
			ReferenceSequenceFile reference,
			int kmer) {
		this.dictionary = dictionary;
		this.reference = reference;
		this.graphf = new DeBruijnReadGraph(kmer, BreakpointDirection.Forward);
		this.graphb = new DeBruijnReadGraph(kmer, BreakpointDirection.Backward);
	}
	private DirectedBreakpointAssembly createDirectedEvidence(AnomolousReadAssembly ara, BreakpointDirection direction) {
		// AnomolousReadAssembly is a placeholder soft-clipped read containing the assembly information
		SoftClipEvidence assembledSC = new SoftClipEvidence(direction, ara);
		return DirectedBreakpointAssembly.create(dictionary, reference, ASSEMBLER_NAME, currentReferenceIndex, currentPosition, BreakpointDirection.Forward,
				assembledSC.getBreakpointSequence(), assembledSC.getBreakpointQuality(),
				ara.getReadBases(), ara.getBaseQualities(),
				ara.getReadCount(), assembledSC.getSoftClipLength());
	}
	@Override
	public Iterable<DirectedBreakpointAssembly> addEvidence(DirectedEvidence evidence) {
		if (evidence == null) return Collections.emptyList();
		BreakpointLocation location = evidence.getBreakpointLocation();
		List<DirectedBreakpointAssembly> result = processUpToExcluding(location.referenceIndex, location.start);
		// add evidence
		if (location.direction == BreakpointDirection.Forward) {
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
	private static void addToGraph(DeBruijnReadGraph graph, PriorityQueue<DirectedEvidence> active, DirectedEvidence evidence) {
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair nrrp = (NonReferenceReadPair)evidence;
			graph.addRead(nrrp.getNonReferenceRead(), false);
			active.add(evidence);
		} else if (evidence instanceof SoftClipEvidence) {
			SoftClipEvidence sc = (SoftClipEvidence)evidence;
			graph.addRead(sc.getSAMRecord(), true);
			active.add(evidence);
		}
	}
	private static void flushUpToExcluding(DeBruijnReadGraph graph, PriorityQueue<DirectedEvidence> active, int referenceIndex, int position) {
		while (active.peek().getBreakpointLocation().referenceIndex < referenceIndex || 
				(active.peek().getBreakpointLocation().referenceIndex == referenceIndex && active.peek().getBreakpointLocation().end < position)) {
			DirectedEvidence evidence = active.poll();
			boolean anchored;
			SAMRecord read;
			if (evidence instanceof NonReferenceReadPair) {
				anchored = false;
				read = ((NonReferenceReadPair)evidence).getNonReferenceRead();
			} else if (evidence instanceof SoftClipEvidence) {
				anchored = true;
				read = ((SoftClipEvidence)evidence).getSAMRecord();
			} else {
				throw new IllegalStateException(String.format("Sanity check failure: unhandled evidence of type %s present in de bruijn graph", evidence.getClass()));
			}
			graph.removeRead(read, anchored);
		}
	}
	@Override
	public Iterable<DirectedBreakpointAssembly> endOfEvidence() {
		return processUpToExcluding(currentReferenceIndex + 1, 0);
	}
	/**
	 * Processes all positions before the given position
	 * @param referenceIndex reference index of position to process until
	 * @param position position to process until
	 * @return add assemblies generated up to the given position
	 */
	private List<DirectedBreakpointAssembly> processUpToExcluding(int referenceIndex, int position) {
		List<DirectedBreakpointAssembly> result = null;
		// keep going till we flush both de bruijn graphs
		// or we reach the stop position
		while (!activef.isEmpty() && !activeb.isEmpty() &&
			!(currentReferenceIndex == referenceIndex && currentPosition == position)) {
			AnomolousReadAssembly araf = graphf.assembleVariant();
			if (araf != null) {
				DirectedBreakpointAssembly evidence = createDirectedEvidence(araf, BreakpointDirection.Forward);
				if (result == null) result = Lists.newArrayList(evidence);
				else result.add(evidence);
			}
			AnomolousReadAssembly arab = graphb.assembleVariant();
			if (arab != null) {
				DirectedBreakpointAssembly evidence = createDirectedEvidence(arab, BreakpointDirection.Backward);
				if (result == null) result = Lists.newArrayList(evidence);
				else result.add(evidence);
			}
			flushUpToExcluding(graphf, activef, currentReferenceIndex, currentPosition);
			flushUpToExcluding(graphb, activeb, currentReferenceIndex, currentPosition);
			currentPosition++;
		}
		currentReferenceIndex = referenceIndex;
		currentPosition = position;
		return result;
	}
}
