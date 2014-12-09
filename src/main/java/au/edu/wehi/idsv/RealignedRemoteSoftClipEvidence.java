package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.ArrayUtils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Soft clip split read in which the realigned soft clip maps to a local coordinate
 * 
 * @author Daniel Cameron
 *
 */
public class RealignedRemoteSoftClipEvidence extends RealignedSoftClipEvidence implements RemoteEvidence {
	private final BreakpointSummary location;
	protected RealignedRemoteSoftClipEvidence(ProcessingContext processContext,
			SAMEvidenceSource source, BreakendDirection direction,
			SAMRecord record, SAMRecord realigned)
			throws CloneNotSupportedException {
		super(processContext, source, direction, record, realigned);
		this.location = super.getBreakendSummary().remoteBreakpoint();
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return location;
	}
	public String getUntemplatedSequence() {
		String seq = super.getUntemplatedSequence();
		if (getRealignedSAMRecord().getReadNegativeStrandFlag()) {
			seq = SequenceUtil.reverseComplement(seq);
		}
		return seq;
	}
	@Override
	public byte[] getBreakendSequence() {
		if (!isBreakendExact()) return null;
		// breakend sequence from this side is untemplated + anchor
		byte[] untemplated = getUntemplatedSequence().getBytes(StandardCharsets.US_ASCII);
		byte[] anchor = getSAMRecord().getReadBases();
		if (super.getBreakendSummary().direction == BreakendDirection.Forward) {
			anchor = ArrayUtils.subarray(anchor, 0, super.getLocalBaseLength());
		} else {
			anchor = ArrayUtils.subarray(anchor, super.getSoftClipLength(), super.getSAMRecord().getReadLength());
		}
		if (getRealignedSAMRecord().getReadNegativeStrandFlag()) {
			SequenceUtil.reverseComplement(anchor);
		}
		if (location.direction == BreakendDirection.Forward) {
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
