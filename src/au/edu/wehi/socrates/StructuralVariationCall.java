package au.edu.wehi.socrates;

import htsjdk.variant.variantcontext.VariantContext;

public class StructuralVariationCall extends VariantContextDirectedBreakpoint {
	protected StructuralVariationCall(ProcessingContext processContext, VariantContext context) {
		super(processContext, context);
	}
}
