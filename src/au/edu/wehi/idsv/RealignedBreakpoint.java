package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

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
		// adjust bounds for microhomologies (but don't overrun contig bounds)
		summary = new BreakpointSummary(
				local.referenceIndex, local.direction,
				local.direction == BreakendDirection.Backward ? local.start : Math.max(local.start - microhomologyLength, 1),
				local.direction == BreakendDirection.Forward ? local.end: Math.min(local.end + microhomologyLength, context.getDictionary().getSequence(local.referenceIndex).getSequenceLength()),
				realigned.getReferenceIndex(), remoteDirection,
				remoteDirection == BreakendDirection.Forward ? remotePosition : Math.max(remotePosition - microhomologyLength, 1),
				remoteDirection == BreakendDirection.Backward ? remotePosition: Math.min(remotePosition + microhomologyLength, context.getDictionary().getSequence(realigned.getReferenceIndex()).getSequenceLength()));
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
		int windowStart = realignDirection == BreakendDirection.Forward ? realignPosition + 1 : realignPosition - 1 - anchoredSequence.length + 1;
		int windowEnd = realignDirection == BreakendDirection.Forward ? realignPosition + 1 + anchoredSequence.length - 1 : realignPosition - 1;
		byte[] referenceBases = getPaddedSubsequenceAt(context, realignReferenceIndex, windowStart, windowEnd);
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
	/**
	 * Gets the subsequence at the given position, padding with Ns if the contig bounds are overrun
	 * @param context
	 * @param realignReferenceIndex
	 * @param windowStart
	 * @param windowEnd
	 * @return
	 */
	private static byte[] getPaddedSubsequenceAt(ProcessingContext context, int realignReferenceIndex, int windowStart, int windowEnd) {
		SAMSequenceRecord contig = context.getDictionary().getSequence(realignReferenceIndex);
		byte[] bases = context.getReference().getSubsequenceAt(contig.getSequenceName(), Math.max(1, windowStart), Math.min(windowEnd, contig.getSequenceLength())).getBases();
		if (windowStart < 1) {
			int toPad = 1 - windowStart;
			byte[] newBases = new byte[bases.length + toPad];
			Arrays.fill(newBases, (byte)'N');
			System.arraycopy(bases, 0, newBases, toPad, bases.length);
			bases = newBases;
		}
		if (windowEnd > contig.getSequenceLength()) {
			int toPad = windowEnd - contig.getSequenceLength();
			byte[] newBases = new byte[bases.length + toPad];
			Arrays.fill(newBases, (byte)'N');
			System.arraycopy(bases, 0, newBases, 0, bases.length);
			bases = newBases;
		}
		return bases;
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
