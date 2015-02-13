package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class RealignedBreakpoint {
	private final BreakpointSummary summary;
	private final String insertedSequence;
	private final boolean exact;
	private RealignedBreakpoint(BreakpointSummary summary, String insertedSequence, boolean exact) {
		this.summary = summary;
		this.insertedSequence = insertedSequence;
		this.exact = exact;
	}
	public static RealignedBreakpoint create(ReferenceSequenceFile reference, SAMSequenceDictionary dictionary, BreakendSummary local, String anchoredSequence, SAMRecord realigned) {
		return create(reference,dictionary, local, anchoredSequence.getBytes(StandardCharsets.US_ASCII), realigned);
	}
	public static RealignedBreakpoint create(ReferenceSequenceFile reference, SAMSequenceDictionary dictionary, BreakendSummary local, byte[] anchoredSequence, SAMRecord realigned) {
		if (realigned.getReadUnmappedFlag()) throw new IllegalArgumentException("Realigned read is not mapped");
		boolean exact = anchoredSequence != null && anchoredSequence.length > 0;
		int ci = local.end - local.start;
		if (exact && ci != 0) {
			throw new IllegalArgumentException("Breakend of " + local.toString() + " not possible for exact breakpoint");
		}
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
		String insertedSequence = insSeq;
		BreakpointSummary summary = new BreakpointSummary(
				local.referenceIndex, local.direction,
				local.start,
				local.end,
				realigned.getReferenceIndex(), remoteDirection,
				remoteDirection == BreakendDirection.Forward ? remotePosition : remotePosition - ci,
				remoteDirection == BreakendDirection.Backward ? remotePosition : remotePosition + ci);
		
		int microhomologyLengthLocalSequenceRemotePosition = 0;
		int microhomologyLengthRemoteSequenceLocalPosition = 0;
		if (insertedSequence.length() == 0 && exact) {
			// only consider microhomology if there is no inserted sequence
			microhomologyLengthLocalSequenceRemotePosition = calculateMicrohomologyLength(reference, dictionary, local.direction, anchoredSequence, realigned.getReferenceIndex(), remoteDirection, remotePosition);
			microhomologyLengthRemoteSequenceLocalPosition = calculateMicrohomologyLength(reference, dictionary, remoteDirection, realigned.getReadBases(), local.referenceIndex, local.direction, local.start);
		}
		if (microhomologyLengthLocalSequenceRemotePosition > 0 || microhomologyLengthRemoteSequenceLocalPosition > 0) {
			summary = new BreakpointSummary(
					summary.referenceIndex, summary.direction,
					summary.direction == BreakendDirection.Backward ? summary.start - microhomologyLengthRemoteSequenceLocalPosition : summary.start - microhomologyLengthLocalSequenceRemotePosition,
					summary.direction == BreakendDirection.Forward ? summary.end + microhomologyLengthRemoteSequenceLocalPosition : summary.end + microhomologyLengthLocalSequenceRemotePosition,
					summary.referenceIndex2, summary.direction2,
					summary.direction2 == BreakendDirection.Forward ? summary.start2 - microhomologyLengthRemoteSequenceLocalPosition : summary.start2 - microhomologyLengthLocalSequenceRemotePosition,
					summary.direction2 == BreakendDirection.Backward ? summary.end2 + microhomologyLengthRemoteSequenceLocalPosition : summary.end2 + microhomologyLengthLocalSequenceRemotePosition) ;
		}
		if (local instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)local;
			if (BreakpointSummary.ByStartStart2EndEnd2.compare(bp, summary) != 0) {
				throw new IllegalArgumentException(String.format("realignment to %s does not match existing breakpoint %s", bp, summary));
			}
		}
		// make sure we don't overrun contig bounds
		summary = new BreakpointSummary(
				summary.referenceIndex, summary.direction,
				ensureOnContig(dictionary, summary.referenceIndex, summary.start),
				ensureOnContig(dictionary, summary.referenceIndex, summary.end),
				summary.referenceIndex2, summary.direction2,
				ensureOnContig(dictionary, summary.referenceIndex2, summary.start2),
				ensureOnContig(dictionary, summary.referenceIndex2, summary.end2));
		return new RealignedBreakpoint(summary, insertedSequence, exact);
	}
	private static int ensureOnContig(SAMSequenceDictionary dict, int referenceIndex, int position) {
		int contigLength = dict.getSequence(referenceIndex).getSequenceLength();
		int adjustedPosition = Math.min(Math.max(position, 1), contigLength); 
		return adjustedPosition;
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
			ReferenceSequenceFile reference,
			SAMSequenceDictionary dictionary,
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
		byte[] referenceBases = getPaddedSubsequenceAt(reference, dictionary, realignReferenceIndex, windowStart, windowEnd);
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
	private static byte[] getPaddedSubsequenceAt(ReferenceSequenceFile reference, SAMSequenceDictionary dictionary, int realignReferenceIndex, int windowStart, int windowEnd) {
		SAMSequenceRecord contig = dictionary.getSequence(realignReferenceIndex);
		byte[] bases = reference.getSubsequenceAt(contig.getSequenceName(), Math.max(1, windowStart), Math.min(windowEnd, contig.getSequenceLength())).getBases();
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
		int c1 = Character.toUpperCase(anchor);
		int c2 = Character.toUpperCase(realign);
		// require exact base matching (any case)
		return c1 == c2 && c1 != 'N';
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
		return isExact() ? summary.end - summary.start : 0;
	}
	public boolean isExact() {
		return exact;
	}
}
