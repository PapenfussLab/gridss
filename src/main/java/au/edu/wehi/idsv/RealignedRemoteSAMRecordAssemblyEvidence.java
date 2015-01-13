package au.edu.wehi.idsv;

import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.ArrayUtils;

public class RealignedRemoteSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence implements RemoteEvidence {
	private final BreakpointSummary bs;
	private final RealignedSAMRecordAssemblyEvidence local;
	public RealignedRemoteSAMRecordAssemblyEvidence(RealignedSAMRecordAssemblyEvidence assembly) {
		super(assembly.getEvidenceSource(), assembly.getSAMRecord(), assembly.getRemoteSAMRecord());
		this.bs = super.getBreakendSummary().remoteBreakpoint();
		this.local = assembly;
	}
	public RealignedSAMRecordAssemblyEvidence asLocal() {
		return local;
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
	public String getEvidenceID() {
		return "R" + super.getEvidenceID();
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
	@Override
	public float getBreakpointQual() {
		return local.getBreakpointQual();
	}
	@Override
	public float getBreakendQual() {
		return local.getBreakendQual();
	}
}
