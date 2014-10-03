package au.edu.wehi.idsv;


/**
 * Annotates breakends with reference allele coverage information
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialCoverageAnnotator implements BreakendAnnotator {
	private final ProcessingContext context;
	private final SequentialReferenceCoverageLookup referenceNormal;
	private final SequentialReferenceCoverageLookup referenceTumour;
	public SequentialCoverageAnnotator(
			ProcessingContext context,
			SequentialReferenceCoverageLookup referenceNormal,
			SequentialReferenceCoverageLookup referenceTumour) {
		this.context = context;
		this.referenceNormal = referenceNormal;
		this.referenceTumour = referenceTumour;
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.BreakendAnnotator#annotate(au.edu.wehi.idsv.VariantContextDirectedEvidence)
	 */
	@Override
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
}
