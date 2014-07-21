package au.edu.wehi.idsv;

import java.util.HashSet;
import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextDirectedBreakpoint extends VariantContextDirectedEvidence implements DirectedBreakpoint {
	public VariantContextDirectedBreakpoint(ProcessingContext processContext, Set<EvidenceSource> sourceSet, VariantContext context, byte[] breakendBaseQual) {
		super(processContext, sourceSet, context, breakendBaseQual);
		assert(super.getBreakendSummary() instanceof BreakpointSummary);
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return (BreakpointSummary)super.getBreakendSummary();
	}
}
