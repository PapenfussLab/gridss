package au.edu.wehi.idsv.alignment;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

/**
 * Determines the length of any inexact breakpoint homology
 * @author Daniel Cameron
 *
 */
public class BreakpointHomology {
	private final int localHomologyLength;
	private final int remoteHomologyLength;
	public BreakpointHomology(int local, int remote) {
		this.localHomologyLength = local;
		this.remoteHomologyLength = remote;
	}
	@Override
	public String toString() {
		return String.format("[%d, %d]", localHomologyLength, remoteHomologyLength);
	}
	/**
	 * Calculates the sequence homology length at the given breakpoint position 
	 * @param lookup reference genome
	 * @param bs breakpoint
	 * @param maxBreakendLength maximum homology length to report
	 * @param additional reference bases to include to accomodate indels
	 * @return breakpoint homology length
	 */
	public static BreakpointHomology calculate(ReferenceLookup lookup, BreakpointSummary bs, String insertedSequence, int maxBreakendLength, int margin) {
		if (bs.start - bs.end != 0 || bs.start2 - bs.end2 != 0) throw new IllegalArgumentException("Breakpoint position must be exact");
		if (insertedSequence == null) insertedSequence = "";
		int seqLength = maxBreakendLength;
		int refLength = maxBreakendLength + insertedSequence.length() + margin;
		if (bs.getEventSize() != null) {
			seqLength = Math.min(seqLength, bs.getEventSize());
			refLength = Math.min(refLength, bs.getEventSize());
		}
		// local           remote
		// ACGTACGT        CCTTAAGG
		//    >                <
		// >>>>                >>>>
		// localSeq           remoteSeq
		//      >>>>       >>>>
		//      localRef   remoteRef
		String localSeq = getAnchorSeq(lookup, bs, refLength);
		String localBsSeq = getAnchorSeq(lookup, bs, seqLength);
		String localRef = getAnchorSeq(lookup, advance(bs, refLength), refLength);
		String remoteSeq = SequenceUtil.reverseComplement(getAnchorSeq(lookup, bs.remoteBreakend(), refLength));
		String remoteBsSeq = SequenceUtil.reverseComplement(getAnchorSeq(lookup, bs.remoteBreakend(), seqLength));
		String remoteRef = SequenceUtil.reverseComplement(getAnchorSeq(lookup, advance(bs.remoteBreakend(), refLength), refLength));
		byte[] breakend = (localBsSeq + insertedSequence + remoteBsSeq).getBytes(StandardCharsets.US_ASCII);
		byte[] local = (localSeq + localRef).getBytes(StandardCharsets.US_ASCII);
		byte[] remote = (remoteRef + remoteSeq).getBytes(StandardCharsets.US_ASCII);
		Aligner aligner = AlignerFactory.create();
		int localHomologyBaseCount = 0;
		int remoteHomologyBaseCount = 0;
		if (breakend != null && breakend.length > 0) {
			if (local != null && local.length > 0) {
				Alignment localAlignment = aligner.align_smith_waterman(breakend, local);
				remoteHomologyBaseCount = remoteBsSeq.length() - SAMRecordUtil.getEndSoftClipLength(TextCigarCodec.decode(localAlignment.getCigar()).getCigarElements());
			}
			if (remote != null && remote.length > 0) {
				Alignment remoteAlignment = aligner.align_smith_waterman(breakend, remote);
				localHomologyBaseCount = localBsSeq.length() - SAMRecordUtil.getStartSoftClipLength(TextCigarCodec.decode(remoteAlignment.getCigar()).getCigarElements());
			}
		}
		return new BreakpointHomology(localHomologyBaseCount, remoteHomologyBaseCount);
	}
	/**
	 * Moves the given breakend forward by the given amount. 
	 */
	private static BreakendSummary advance(BreakendSummary bs, int bases) {
		int offset = bases;
		if (bs.direction == BreakendDirection.Backward) {
			offset *= -1;
		}
		return new BreakendSummary(bs.referenceIndex, bs.direction, bs.start + offset, bs.end + offset);
	}
	private static String getAnchorSeq(final ReferenceLookup lookup, final BreakendSummary bs, final int length) {
		final SAMSequenceRecord refseq = lookup.getSequenceDictionary().getSequence(bs.referenceIndex);
		int start;
		int end;
		if (bs.direction == BreakendDirection.Forward) {
			end = bs.start;
			start = end - length + 1;
		} else {
			start = bs.start;
			end = start + length - 1;
		}
		start = Math.max(1, start);
		end = Math.min(refseq.getSequenceLength(), end);
		byte[] bseq = lookup.getSubsequenceAt(refseq.getSequenceName(), start, end).getBases();
		if (bs.direction == BreakendDirection.Backward) {
			SequenceUtil.reverseComplement(bseq);
		}
		String seq = new String(bseq);
		return seq;
	}
	public int getLocalHomologyLength() {
		return localHomologyLength;
	}
	public int getRemoteHomologyLength() {
		return remoteHomologyLength;
	}
}
