package au.edu.wehi.socrates;

import org.broadinstitute.variant.variantcontext.VariantContext;

public class StructuralVariationCall extends VariantContextDirectedBreakpoint {
	protected StructuralVariationCall(ProcessingContext processContext, VariantContext context) {
		super(processContext, context);
	}
}
