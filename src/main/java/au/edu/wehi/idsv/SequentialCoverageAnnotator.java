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
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(context, variant);
		annotateReferenceCounts(builder, loc.referenceIndex, loc.start + (loc.direction == BreakendDirection.Forward ? 0 : 1));
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
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		if (!it.hasNext()) return endOfData();
		return annotate(it.next());
	}
}
