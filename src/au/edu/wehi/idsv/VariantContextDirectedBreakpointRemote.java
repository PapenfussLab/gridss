package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextDirectedBreakpointRemote extends VariantContextDirectedBreakpoint implements RemoteEvidence {
	private BreakpointSummary location;
	public VariantContextDirectedBreakpointRemote(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(processContext, source, context);
		this.location = super.getBreakendSummary().remoteBreakpoint();
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return location;
	}
}
