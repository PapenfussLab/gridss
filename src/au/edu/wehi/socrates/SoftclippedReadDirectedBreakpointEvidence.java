package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.variant.variantcontext.VariantContext;

import au.edu.wehi.socrates.util.SAMRecordSummary;

public class SoftclippedReadDirectedBreakpointEvidence implements DirectedBreakpoint {
	private BreakpointDirection direction;
	private SAMRecord record;
	private SAMRecordSummary summary;
	public SoftclippedReadDirectedBreakpointEvidence(BreakpointDirection direction, SAMRecord record) {
		this.direction = direction;
		this.record = record;
		this.summary = new SAMRecordSummary(record);
	}
	@Override
	public String getBreakpointID() {
		return String.format("r_%d_%s", record.getFirstOfPairFlag(), record.getReadName());
	}
	@Override
	public BreakpointDirection getBreakpointDirection() {
		return direction;
	}
	@Override
	public byte[] getBreakpointSequence() {
		return direction == BreakpointDirection.Forward ? summary.getTailClipSequence() : summary.getHeadClipSequence();
	}
	@Override
	public byte[] getBreakpointQuality() {
		return direction == BreakpointDirection.Forward ? summary.getTailClipQuality() : summary.getHeadClipQuality();
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public int getSoftClipLength() {
		return getBreakpointSequence().length;
	}
	public float getAlignedPercentIdentity() {
		return summary.getAlignedPercentIdentity();
	}
	public float getAverageClipQuality() {
		return direction == BreakpointDirection.Forward ? summary.getAvgTailClipQuality() : summary.getAvgHeadClipQuality();
	}
	@Override
	public VariantContext getVariantContext() {
		return null;
	}
}
