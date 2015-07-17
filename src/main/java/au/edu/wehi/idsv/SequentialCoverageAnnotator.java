package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.List;


/**
 * Annotates breakends with reference allele coverage information
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialCoverageAnnotator implements Iterator<VariantContextDirectedEvidence>, BreakendAnnotator {
	private final ProcessingContext context;
	private final List<ReferenceCoverageLookup> reference;
	private final Iterator<VariantContextDirectedEvidence> it;
	public SequentialCoverageAnnotator(
			ProcessingContext context,
			Iterator<VariantContextDirectedEvidence> it,
			List<ReferenceCoverageLookup> reference) {
		this.it = it;
		this.context = context;
		this.reference = reference;
	}
	public VariantContextDirectedEvidence annotate(VariantContextDirectedEvidence variant) {
		BreakendSummary loc = variant.getBreakendSummary();
		int referenceIndex = loc.referenceIndex;
		// TODO: start or end based on call position
		int positionImmediatelyBeforeBreakend = loc.start - (loc.direction == BreakendDirection.Forward ? 0 : 1);
		int[] reads = new int[reference.size()];
		int[] spans = new int[reference.size()];
		for (int i = 0; i < reference.size(); i++) {
			ReferenceCoverageLookup r = reference.get(i);
			if (r != null) {
				reads[i] = r.readsSupportingNoBreakendAfter(referenceIndex, positionImmediatelyBeforeBreakend);
				spans[i] = r.readPairsSupportingNoBreakendAfter(referenceIndex, positionImmediatelyBeforeBreakend);
			}
		}
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context, variant);
		builder.referenceReads(reads);
		builder.referenceSpanningPairs(spans);
		return (VariantContextDirectedEvidence)builder.make();
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}
	@Override
	public VariantContextDirectedEvidence next() {
		return annotate(it.next());
	}
}
