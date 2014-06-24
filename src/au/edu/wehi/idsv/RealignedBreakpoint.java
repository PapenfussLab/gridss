package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import au.edu.wehi.idsv.vcf.VcfAttributes;

public class RealignedBreakpoint {
	private final BreakpointSummary summary;
	private final String insertedSequence;
	public RealignedBreakpoint(ProcessingContext context, BreakendSummary local, String anchoredSequence, SAMRecord realigned) {
		this(context, local, anchoredSequence.getBytes(StandardCharsets.US_ASCII), realigned);
	}
	public RealignedBreakpoint(ProcessingContext context, BreakendSummary local, byte[] anchoredSequence, SAMRecord realigned) {
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
		int microhomologyLength = 0;
		if (insertedSequence.length() == 0) {
			// only consider microhomology if there is no inserted sequence
			microhomologyLength = calculateMicrohomologyLength(context, local.direction, anchoredSequence, realigned.getReferenceIndex(), remoteDirection, remotePosition);
		}
		// adjust bounds for microhomologies
		summary = new BreakpointSummary(
				local.referenceIndex, local.direction,
				local.direction == BreakendDirection.Backward ? local.start : local.start - microhomologyLength,
				local.direction == BreakendDirection.Forward ? local.end: local.end + microhomologyLength,
				realigned.getReferenceIndex(), remoteDirection,
				remoteDirection == BreakendDirection.Forward ? remotePosition : remotePosition - microhomologyLength,
				remoteDirection == BreakendDirection.Backward ? remotePosition: remotePosition + microhomologyLength,
				new EvidenceMetrics());
		EvidenceMetrics evidence = summary.evidence;
		evidence.add(local.evidence);
		evidence.set(VcfAttributes.REALIGN_MAX_LENGTH, realigned.getReadLength());
		evidence.set(VcfAttributes.REALIGN_TOTAL_LENGTH, realigned.getReadLength());
		//evidence.set(EvidenceAttributes.REALIGN_MAX_MAPQ, realigned.getMappingQuality());
		evidence.set(VcfAttributes.REALIGN_TOTAL_MAPQ, realigned.getMappingQuality());
		if (local instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)local;
			if (BreakpointSummary.ByStartStart2EndEnd2.compare(bp, summary) != 0) {
				throw new IllegalArgumentException(String.format("realignment to %s does not match existing breakpoint %s", bp, summary));
			}
		}
	}
	/**
	 * Determines how many of the anchor bases match the reference at the realignment location
	 * @param context
	 * @param anchorDirection direction of anchor breakend
	 * @param anchoredSequence anchored sequence
	 * @param remoteLocation realignment location
	 * @param baseOffset
	 * @return
	 */
	private static int calculateMicrohomologyLength(
			ProcessingContext context,
			BreakendDirection anchorDirection,
			byte[] anchoredSequence,
			Integer realignReferenceIndex,
			BreakendDirection realignDirection,
			int realignPosition) {
		// Micro-homology detection:
		// anchored sequence: TTTAAAAA>
		// realign reference: GGGAAAAA
		if (realignReferenceIndex == null || realignReferenceIndex < 0) return 0;
		byte[] referenceBases = context.getReference()
			.getSubsequenceAt(context.getDictionary().getSequence(realignReferenceIndex).getSequenceName(),
				realignDirection == BreakendDirection.Forward ? realignPosition + 1 : realignPosition - 1 - anchoredSequence.length + 1,
				realignDirection == BreakendDirection.Forward ? realignPosition + 1 + anchoredSequence.length - 1 : realignPosition - 1)
			.getBases();
		// set up traversal iterators
		int firstAnchorBase = 0;
		int anchorInc = 1;
		if (anchorDirection == BreakendDirection.Forward) {
			firstAnchorBase = anchoredSequence.length - 1;
			anchorInc = -1;
		}
		int firstRealignBase = 0;
		int realignInc = 1;
		if (realignDirection == BreakendDirection.Backward) {
			firstRealignBase = anchoredSequence.length - 1;
			realignInc = -1;
		}
		int i = 0;
		while (i < anchoredSequence.length && basesMatch(anchoredSequence[firstAnchorBase + anchorInc * i], referenceBases[firstRealignBase + realignInc * i])) {
			i++;
		}
		return i;
	}
	private static boolean basesMatch(byte anchor, byte realign) {
		// require exact base matching (any case)
		return Character.toUpperCase(anchor) == Character.toUpperCase(realign);
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
	public int getMicroHomologyLength() {
		return summary.end - summary.start;
	}
}
