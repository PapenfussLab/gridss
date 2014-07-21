package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Iterator;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.vcf.PassFiltersIterator;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Iterates through evidence for a breakpoint at 
 * @author Daniel Cameron
 *
 */
public class DirectedEvidenceIterator extends AbstractIterator<DirectedEvidence> {
	private static final int INITIAL_BUFFER_SIZE = 32;
	/**
	 * Number of fragment lengths to advance before matching breakend evidence with the realignment SAMRecord.
	 * This delay is required as the realignments are ordered by breakend position, whereas SAMRecord evidence
	 * is traversed according to the SAMRecord start coordinate.
	 */
	private static final float WINDOW_LAG_FRAG_LENGTHS = 1;
	private final ProcessingContext processContext;
	private final EvidenceSource source;
	private final PeekingIterator<SAMRecord> svIterator;
	private final PeekingIterator<VariantContext> vcfIterator;
	private final SequentialNonReferenceReadPairFactory mateFactory;
	private final SequentialRealignedBreakpointFactory realignFactory;
	private final PriorityQueue<DirectedEvidence> calls = new PriorityQueue<DirectedEvidence>(INITIAL_BUFFER_SIZE, DirectedEvidenceOrder.ByNatural);
	private long currentLinearPosition = -1;
	/**
	 * Iterates over breakpoint evidence
	 * @param sv structural variation supporting reads. Must be sorted by coordinate
	 * @param mate mate reads of structural variation supporting reads. Must be sorted by mate coordinate
	 * @param realign breakpoint realignments. Must be sorted by source @see DirectedBreakpoint start coordinate 
	 * @param vcf Assembly-based @see DirectedBreakpoint. Must be sorted by coordinate
	 */
	public DirectedEvidenceIterator(
			ProcessingContext processContext,
			EvidenceSource source,
			Iterator<SAMRecord> sv,
			Iterator<SAMRecord> mate,
			Iterator<SAMRecord> realign,
			Iterator<VariantContext> vcf) {
		this.processContext = processContext;
		this.source = source;
		this.svIterator = sv == null ? null : Iterators.peekingIterator(sv);
		this.vcfIterator = vcf == null ? null : Iterators.peekingIterator(new PassFiltersIterator<VariantContext>(vcf));
		this.mateFactory = mate == null ? null : new SequentialNonReferenceReadPairFactory(Iterators.peekingIterator(mate));
		this.realignFactory = realign == null ? null : new SequentialRealignedBreakpointFactory(Iterators.peekingIterator(realign));
	}
	@Override
	protected DirectedEvidence computeNext() {
		do {
			if (!calls.isEmpty() && processContext.getLinear().getStartLinearCoordinate(calls.peek().getBreakendSummary()) < currentLinearPosition - WINDOW_LAG_FRAG_LENGTHS * source.getMetrics().getMaxFragmentSize()) {
				return fixCall(calls.poll());
			}
		} while (advance());
		// no more input: flush our buffer
		while (!calls.isEmpty()) return fixCall(calls.poll());
		return endOfData();
	}
	/**
	 * Need to match realigned record at output time since input BAM processing is in order of alignment start,
	 * not putative breakpoint position
	 * @param evidence call to match realignment record for
	 * @return call including realignment mapping evidence if applicable
	 */
	private DirectedEvidence fixCall(DirectedEvidence evidence) {
		if (realignFactory == null) return evidence;
		if (evidence instanceof SoftClipEvidence) {
			evidence = SoftClipEvidence.create((SoftClipEvidence)evidence, realignFactory.findRealignedSAMRecord(evidence));
		} else if (evidence instanceof VariantContextDirectedEvidence) {
			VariantContextDirectedEvidence assembly = (VariantContextDirectedEvidence)evidence;
			SAMRecord realigned = realignFactory.findRealignedSAMRecord(assembly);
			evidence = AssemblyBuilder.incorporateRealignment(processContext, assembly, realigned);
		}
		return evidence;
	}
	/**
	 * Advances the input streams to the next
	 * @return true if records were consumed, false if no more input exists
	 */
	private boolean advance() {
		long nextPos = Long.MAX_VALUE;
		if (svIterator != null && svIterator.hasNext()) {
			nextPos = Math.min(nextPos, processContext.getLinear().getLinearCoordinate(svIterator.peek().getReferenceIndex(), svIterator.peek().getAlignmentStart()));
		}
		if (vcfIterator != null && vcfIterator.hasNext()) {
			nextPos = Math.min(nextPos, processContext.getLinear().getLinearCoordinate(vcfIterator.peek().getChr(), vcfIterator.peek().getStart()));
		}
		if (nextPos == Long.MAX_VALUE) return false; // no records
		while (svIterator != null && svIterator.hasNext() && processContext.getLinear().getLinearCoordinate(svIterator.peek().getReferenceIndex(), svIterator.peek().getAlignmentStart()) == nextPos) {
			processRead(svIterator.next());
		}
		while (vcfIterator != null && vcfIterator.hasNext() && processContext.getLinear().getLinearCoordinate(vcfIterator.peek().getChr(), vcfIterator.peek().getStart()) == nextPos) {
			processVariant(vcfIterator.next());
		}
		currentLinearPosition = nextPos;
		return true;
	}
	private void processVariant(VariantContext variant) {
		IdsvVariantContext managedContext = IdsvVariantContext.create(processContext, source, variant);
		if (managedContext instanceof DirectedEvidence) {
			calls.add((DirectedEvidence)managedContext);
		}
	}
	private void processRead(SAMRecord record) {
		if (SAMRecordUtil.getStartSoftClipLength(record) > 0) {
			SoftClipEvidence sce = new SoftClipEvidence(processContext, (SAMEvidenceSource)source, BreakendDirection.Backward, record);
			if (processContext.getSoftClipParameters().meetsEvidenceCritera(sce)) {
				calls.add(sce);
			}
		}
		if (SAMRecordUtil.getEndSoftClipLength(record) > 0) {
			SoftClipEvidence sce = new SoftClipEvidence(processContext, (SAMEvidenceSource)source, BreakendDirection.Forward, record);
			if (processContext.getSoftClipParameters().meetsEvidenceCritera(sce)) {
				calls.add(sce);
			}
		}
		if (mateFactory != null && SAMRecordUtil.isPartOfNonReferenceReadPair(record)) {
			NonReferenceReadPair pair = mateFactory.createNonReferenceReadPair(record, (SAMEvidenceSource)source);
			if (pair != null && pair.isValid()) {
				calls.add(pair);
			}
		}
	}
}
