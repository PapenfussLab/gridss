package au.edu.wehi.socrates;

import org.broadinstitute.variant.variantcontext.VariantContext;

public class VariantContextDirectedBreakpoint extends VariantContext implements DirectedBreakpoint {
	protected VariantContextDirectedBreakpoint(VariantContext context) {
		super(context);
		throw new RuntimeException("NYI");
	}
	@Override
	public BreakpointDirection getBreakpointDirection() {
		throw new RuntimeException("NYI");
	}
	@Override
	public byte[] getBreakpointSequence() {
		throw new RuntimeException("NYI");
	}
	@Override
	public VariantContext getVariantContext() {
		throw new RuntimeException("NYI");
	}
	@Override
	public String getBreakpointID() {
		return getID();
	}
	@Override
	public byte[] getBreakpointQuality() {
		// TODO Auto-generated method stub
		return null;
	}
}
