package au.edu.wehi.socrates;

import java.util.PriorityQueue;

import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;

/**
 * Annotates sorted breakends with the evidence supporting the call
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialBreakendAnnotator {
	private final ProcessingContext context;
	private final PeekingIterator<DirectedBreakpoint> evidence;
	private final SequentialReferenceCoverageLookup reference;
	private final PriorityQueue<DirectedEvidence> activeEvidence = new PriorityQueue<DirectedEvidence>(new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return BreakendSummary.ByEndStart.compare(arg0.getBreakendSummary(), arg1.getBreakendSummary());
		}
	});
	private int currentReferenceIndex = -1;
	public SequentialBreakendAnnotator(
			ProcessingContext context,
			SequentialReferenceCoverageLookup reference,
			PeekingIterator<DirectedBreakpoint> evidence) {
		this.context = context;
		this.evidence = evidence;
		this.reference = reference;
	}
	public VariantContextDirectedBreakpoint annotate(VariantContextDirectedBreakpoint variant) {
		BreakendSummary loc = variant.getBreakendSummary();
		if (currentReferenceIndex != loc.referenceIndex) {
			activeEvidence.clear();
			if (currentReferenceIndex > loc.referenceIndex) {
				throw new IllegalArgumentException(String.format("Sanity check failure: variants must be sorted by start position"));
			}
			currentReferenceIndex = loc.referenceIndex;
		}
		flushBefore(loc.referenceIndex, loc.start);
		addUntil(loc.referenceIndex, loc.end);
		
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(context, variant.getBreakendSummary())
			.referenceReads(reference.readsSupportingNoBreakendAfter(loc.referenceIndex, loc.start + (loc.direction == BreakendDirection.Forward ? 0 : 1)))
			.referenceSpanningPairs(reference.readPairsSupportingNoBreakendAfter(loc.referenceIndex, loc.start + (loc.direction == BreakendDirection.Forward ? 0 : 1)));
		for (DirectedEvidence evidence : activeEvidence) {
			builder.evidence(evidence);
		}
		return builder.make();
	}
	private void addUntil(int referenceIndex, int position) {
		while (evidence.hasNext() && evidence.peek().getBreakendSummary().referenceIndex < referenceIndex) evidence.next();
		while (evidence.hasNext() && evidence.peek().getBreakendSummary().referenceIndex == referenceIndex &&
				evidence.peek().getBreakendSummary().start <= position) {
			activeEvidence.add(evidence.next());
		}
	}
	private void flushBefore(int referenceIndex, int start) {
		while (!activeEvidence.isEmpty() && activeEvidence.peek().getBreakendSummary().end < start) activeEvidence.poll();
	}
}
