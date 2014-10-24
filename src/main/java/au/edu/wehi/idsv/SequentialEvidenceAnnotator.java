package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Annotates sorted breakends with the evidence supporting the call
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialEvidenceAnnotator implements BreakendAnnotator {
	private final ProcessingContext context;
	private final PeekingIterator<DirectedEvidence> evidence;
	private final PriorityQueue<DirectedEvidence> activeEvidence = new PriorityQueue<DirectedEvidence>(1024, DirectedEvidence.ByEndStart);
	private int currentReferenceIndex = -1;
	
	public SequentialEvidenceAnnotator(
			ProcessingContext context,
			Iterator<DirectedEvidence> iterator) {
		this.context = context;
		this.evidence = Iterators.peekingIterator(iterator);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.BreakendAnnotator#annotate(au.edu.wehi.idsv.VariantContextDirectedEvidence)
	 */
	@Override
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
		
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(context, variant);
		
		for (DirectedEvidence evidence : activeEvidence) {
			BreakendSummary evidenceLoc = evidence.getBreakendSummary();
			if (evidenceLoc.overlaps(loc)) {
				builder.addEvidence(evidence);
			}
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
