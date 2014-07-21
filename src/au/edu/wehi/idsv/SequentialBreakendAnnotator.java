package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.Iterators;
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
	private final PeekingIterator<DirectedEvidence> evidence;
	private final SequentialReferenceCoverageLookup referenceNormal;
	private final SequentialReferenceCoverageLookup referenceTumour;
	private final PriorityQueue<DirectedEvidence> activeEvidence = new PriorityQueue<DirectedEvidence>(1024, new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return BreakendSummary.ByEndStart.compare(arg0.getBreakendSummary(), arg1.getBreakendSummary());
		}
	});
	private int currentReferenceIndex = -1;
	
	public SequentialBreakendAnnotator(
			ProcessingContext context,
			SequentialReferenceCoverageLookup referenceNormal,
			SequentialReferenceCoverageLookup referenceTumour,
			Iterator<DirectedEvidence> iterator) {
		this.context = context;
		this.evidence = Iterators.peekingIterator(iterator);
		this.referenceNormal = referenceNormal;
		this.referenceTumour = referenceTumour;
	}
	public VariantContextDirectedEvidence annotate(VariantContextDirectedEvidence variant) {
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
		
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(context, variant.getBreakendSummary());
		annotateReferenceCounts(builder, loc.referenceIndex, loc.start + (loc.direction == BreakendDirection.Forward ? 0 : 1));
		
		for (DirectedEvidence evidence : activeEvidence) {
			BreakendSummary evidenceLoc = evidence.getBreakendSummary();
			if (evidenceLoc.overlaps(loc)) {
				builder.addEvidence(evidence);
			}
		}
		return builder.make();
	}
	private void annotateReferenceCounts(StructuralVariationCallBuilder builder, int referenceIndex, int positionImmediatelyBeforeBreakend) {
		int normalReads = 0;
		int normalSpans = 0;
		int tumourReads = 0;
		int tumourSpans = 0;
		if (referenceNormal != null) {
			normalReads = referenceNormal.readsSupportingNoBreakendAfter(referenceIndex, positionImmediatelyBeforeBreakend);
			normalSpans = referenceNormal.readPairsSupportingNoBreakendAfter(referenceIndex, positionImmediatelyBeforeBreakend);
		}
		if (referenceTumour != null) {
			tumourReads = referenceTumour.readsSupportingNoBreakendAfter(referenceIndex, positionImmediatelyBeforeBreakend);
			tumourSpans = referenceTumour.readPairsSupportingNoBreakendAfter(referenceIndex, positionImmediatelyBeforeBreakend);
		}
		builder.referenceReads(normalReads, tumourReads);
		builder.referenceSpanningPairs(normalSpans, tumourSpans);
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
