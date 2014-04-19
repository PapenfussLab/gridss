package au.edu.wehi.socrates;

import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMRecord;

public class SoftclipReadProcessor implements SequentialDirectedBreakpointEvidenceProcessor {
	private List<SoftClipSupport> reads = new ArrayList<SoftClipSupport>(); 
	private final int minMapq;
	private final int minSoftClipLength;
	private final float minAlignedPercentIdentity;
	private final float minAverageSoftClipBaseQuality;
	public SoftclipReadProcessor(int minMapq, int minSoftClipLength, float minAlignedPercentIdentity, float minAverageSoftClipBaseQuality) {
		this.minMapq = minMapq;
		this.minSoftClipLength = minSoftClipLength;
		this.minAlignedPercentIdentity = minAlignedPercentIdentity;
		this.minAverageSoftClipBaseQuality = minAverageSoftClipBaseQuality;
	}
	@Override
	public Iterable<DirectedBreakpoint> getEvidenceAtPosition(int position) {
		// TODO: why does Socrates v1 dedup in BAMStratifier?
		// doesn't MarkDuplicates do this already?
		
		// TODO: should we at least restrict to only unique sequences? 
		List<DirectedBreakpoint> dp = new ArrayList<DirectedBreakpoint>();
		for (SoftClipSupport sss : reads) {
			SoftclippedReadDirectedBreakpointEvidence bp = new SoftclippedReadDirectedBreakpointEvidence(sss.direction, sss.record);
			if (sss.record.getMappingQuality() > minMapq &&
					bp.getSoftClipLength() > minSoftClipLength &&
					bp.getAlignedPercentIdentity() > minAlignedPercentIdentity &&
					bp.getAverageClipQuality() > minAverageSoftClipBaseQuality) {
				dp.add(bp);
			}
		}
		return dp;
	}
	private class SoftClipSupport {
		public SoftClipSupport(BreakpointDirection direction, SAMRecord record) {
			this.direction = direction;
			this.record = record;
		}
		public final BreakpointDirection direction;
		public final SAMRecord record;
	}
	@Override
	public void addSoftclipSupportingForwardBreakpoint(SAMRecord record) {
		reads.add(new SoftClipSupport(BreakpointDirection.Forward, record));
	}
	@Override
	public void addSoftclipSupportingBackwardBreakpoint(SAMRecord record) {
		reads.add(new SoftClipSupport(BreakpointDirection.Backward, record));
	}
	@Override
	public void addReadPairSupportingForwardBreakpoint(NonReferenceReadPair pair) {
	}
	@Override
	public void addReadPairSupportingBackwardBreakpoint(NonReferenceReadPair pair) {
	}
}
