package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.ArrayUtils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class RealignedRemoteSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence implements RemoteEvidence {
	private final BreakpointSummary bs;
	public RealignedRemoteSAMRecordAssemblyEvidence(ProcessingContext processContext, AssemblyEvidenceSource source, SAMRecord assembly, SAMRecord realigned) {
		super(processContext, source, assembly, realigned);
		this.bs = super.getBreakendSummary().remoteBreakpoint();
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return bs;
	}
	public String getUntemplatedSequence() {
		String seq = super.getUntemplatedSequence();
		if (getRemoteSAMRecord().getReadNegativeStrandFlag()) {
			seq = SequenceUtil.reverseComplement(seq);
		}
		return seq;
	}
	@Override
	public byte[] getBreakendSequence() {
		if (!isBreakendExact()) return null;
		// breakend sequence from this side is untemplated + anchor
		byte[] untemplated = getUntemplatedSequence().getBytes(StandardCharsets.US_ASCII);
		byte[] anchor = getAssemblyAnchorSequence();
		if (getRemoteSAMRecord().getReadNegativeStrandFlag()) {
			SequenceUtil.reverseComplement(anchor);
		}
		if (bs.direction == BreakendDirection.Forward) {
			return ArrayUtils.addAll(untemplated, anchor);
		} else {
			return ArrayUtils.addAll(anchor, untemplated);
		}
	}
	@Override
	public String toString() {
		return "R" + super.toString();
	}
	@Override
	public int getLocalBaseLength() {
		return super.getRemoteBaseLength();
	}
	@Override
	public int getLocalMapq() {
		return super.getRemoteMapq();
	}
	@Override
	public int getLocalMaxBaseQual() {
		return super.getRemoteMaxBaseQual();
	}
	@Override
	public int getLocalTotalBaseQual() {
		return super.getRemoteTotalBaseQual();
	}
	@Override
	public int getRemoteBaseCount() {
		return super.getLocalBaseLength();
	}
	@Override
	public int getRemoteBaseLength() {
		return super.getLocalBaseLength();
	}
	@Override
	public int getRemoteMapq() {
		return super.getLocalMapq();
	}
	@Override
	public int getRemoteMaxBaseQual() {
		return super.getLocalMaxBaseQual();
	}
	@Override
	public int getRemoteTotalBaseQual() {
		return super.getLocalTotalBaseQual();
	}
}
