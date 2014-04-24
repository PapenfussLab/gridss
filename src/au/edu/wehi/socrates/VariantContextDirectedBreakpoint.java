package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.variant.variantcontext.VariantContext;

public abstract class VariantContextDirectedBreakpoint extends VariantContext implements RealignedDirectedBreakpoint {
	private final int referenceIndex;
	private SAMRecord realigned;
	protected VariantContextDirectedBreakpoint(VariantContext context, int referenceIndex) {
		super(context);
		this.referenceIndex = referenceIndex;
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
	public String getEvidenceID() {
		return getID();
	}
	@Override
	public byte[] getBreakpointQuality() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public int getReferenceIndex() {
		return referenceIndex;
	}
	@Override
	public long getWindowStart() {
		return getStart();
	}
	@Override
	public long getWindowEnd() {
		return getEnd();
	}
	@Override
	public void setRealigned(SAMRecord realigned) {
		this.realigned = realigned;		
	}
	@Override
	public SAMRecord getRealigned() {
		return realigned;
	}
}
