package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import au.edu.wehi.socrates.vcf.EvidenceAttributes;
import au.edu.wehi.socrates.vcf.VcfConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class RealignedBreakpoint {
	private final BreakpointSummary summary;
	private final String insertedSequence;
	public RealignedBreakpoint(BreakendSummary local, SAMRecord realigned) {
		if (realigned.getReadUnmappedFlag()) throw new IllegalArgumentException("Realigned read is not mapped");
		BreakendDirection remoteDirection;
		int remotePosition;
		int insertedSequenceLength;
		String insSeq;
		if ((local.direction == BreakendDirection.Forward && realigned.getReadNegativeStrandFlag()) ||
				(local.direction == BreakendDirection.Backward && !realigned.getReadNegativeStrandFlag())) {
			// negative strand template match means we flip the expected direction
			// since breakend sequence is on the +ve strand,
			// realigned forward breakends on +ve stand would indicate backward breakend
			// realigned backward breakends on +ve stand would indicate forward breakend
			remoteDirection = BreakendDirection.Forward;
			remotePosition = realigned.getAlignmentEnd();
			insertedSequenceLength = SAMRecordUtil.getEndSoftClipLength(realigned);
			insSeq = new String(Arrays.copyOfRange(realigned.getReadBases(), realigned.getReadLength() - insertedSequenceLength, realigned.getReadLength()), StandardCharsets.US_ASCII);
		} else {
			// ACGT. => CGT breakpoint, -> ACGT.---.CGT
			remoteDirection = BreakendDirection.Backward;
			remotePosition = realigned.getAlignmentStart();
			insertedSequenceLength = SAMRecordUtil.getStartSoftClipLength(realigned);
			insSeq = new String(Arrays.copyOfRange(realigned.getReadBases(), 0, insertedSequenceLength), StandardCharsets.US_ASCII);
		}
		if (realigned.getReadNegativeStrandFlag()) {
			insSeq = SequenceUtil.reverseComplement(insSeq);
		}
		insertedSequence = insSeq;
		BreakendSummary remoteLocation = new BreakendSummary(realigned.getReferenceIndex(), remoteDirection, remotePosition, remotePosition, null);
		summary = new BreakpointSummary(local, remoteLocation, local.evidence == null ? new EvidenceMetrics() : local.evidence);
		EvidenceMetrics evidence = summary.evidence;
		evidence.set(EvidenceAttributes.REALIGN_MAX_LENGTH, realigned.getReadLength());
		evidence.set(EvidenceAttributes.REALIGN_TOTAL_LENGTH, realigned.getReadLength());
		//evidence.set(EvidenceAttributes.REALIGN_MAX_MAPQ, realigned.getMappingQuality());
		evidence.set(EvidenceAttributes.REALIGN_TOTAL_MAPQ, realigned.getMappingQuality());
		if (local instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)local;
			if (BreakpointSummary.ByStartStart2EndEnd2.compare(bp, summary) != 0) {
				throw new IllegalArgumentException(String.format("realignment to %s does not match existing breakpoint %s", bp, summary));
			}
		}
	}
	public BreakpointSummary getBreakpointSummary() {
		return summary;
	}
	public int getInsertedSequenceLength() {
		return insertedSequence.length();
	}
	/**
	 * Returns the untemplated breakpoint sequence not mapped to either the local
	 * or realigned mapping location.
	 * @return
	 */
	public String getInsertedSequence() {
		return insertedSequence;
	}
}
