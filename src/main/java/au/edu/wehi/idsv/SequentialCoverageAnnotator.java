package au.edu.wehi.idsv;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;


/**
 * Annotates breakends with reference allele coverage information
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialCoverageAnnotator extends AbstractIterator<VariantContextDirectedEvidence> implements BreakendAnnotator {
	private final ProcessingContext context;
	private final ReferenceCoverageLookup referenceNormal;
	private final ReferenceCoverageLookup referenceTumour;
	private final Iterator<VariantContextDirectedEvidence> it;
	public SequentialCoverageAnnotator(
			ProcessingContext context,
			Iterator<VariantContextDirectedEvidence> it,
			ReferenceCoverageLookup referenceNormal,
			ReferenceCoverageLookup referenceTumour) {
		this.it = it;
		this.context = context;
		this.referenceNormal = referenceNormal;
		this.referenceTumour = referenceTumour;
	}
	public VariantContextDirectedEvidence annotate(VariantContextDirectedEvidence variant) {
		BreakendSummary loc = variant.getBreakendSummary();
		int referenceIndex = loc.referenceIndex;
		// TODO: start or end based on call position
		int positionImmediatelyBeforeBreakend = loc.start - (loc.direction == BreakendDirection.Forward ? 0 : 1);
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
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context, variant);
		builder.referenceReads(normalReads, tumourReads);
		builder.referenceSpanningPairs(normalSpans, tumourSpans);
		return (VariantContextDirectedEvidence)builder.make();
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		if (!it.hasNext()) return endOfData();
		return annotate(it.next());
	}
}
