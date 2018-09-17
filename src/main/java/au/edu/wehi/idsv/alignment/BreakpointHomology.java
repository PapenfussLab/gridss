package au.edu.wehi.idsv.alignment;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;

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
	 * @param margin additional reference bases to include to accommodate indels
	 * @return breakpoint homology length
	 */
	public static BreakpointHomology calculate(ReferenceLookup lookup, BreakpointSummary bs, String insertedSequence, int maxBreakendLength, int margin) {
		if (bs.start - bs.end != 0 || bs.start2 - bs.end2 != 0) {
			throw new IllegalArgumentException("Breakpoint position must be exact");
		}
		if (insertedSequence == null) {
			insertedSequence = "";
		}
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
		String strBreakend = localBsSeq + insertedSequence + remoteBsSeq;
		String strLocal = localSeq + localRef;
		String strRemote = remoteRef + remoteSeq;
		byte[] breakend = strBreakend.getBytes(StandardCharsets.US_ASCII);
		byte[] local = strLocal.getBytes(StandardCharsets.US_ASCII);
		byte[] remote = strRemote.getBytes(StandardCharsets.US_ASCII);
		Aligner aligner = AlignerFactory.create();
		int localHomologyBaseCount = 0;
		int remoteHomologyBaseCount = 0;
		if (breakend != null && breakend.length > 0) {
			if (local != null && local.length > 0) {
				Alignment localAlignment = aligner.align_smith_waterman(breakend, local);
				List<CigarElement> cigar = TextCigarCodec.decode(localAlignment.getCigar()).getCigarElements();
				// We are defining a homology as the number of bases mapped on the other side
				// inserted sequence means the number of bases consumed can be negative
				remoteHomologyBaseCount = Math.max(0, remoteBsSeq.length() - SAMRecordUtil.getEndSoftClipLength(cigar));
				if (SAMRecordUtil.getStartSoftClipLength(cigar) > 0) {
					// anchor is not aligned - something went wrong
					remoteHomologyBaseCount = 0;
				}
			}
			if (remote != null && remote.length > 0) {
				Alignment remoteAlignment = aligner.align_smith_waterman(breakend, remote);
				List<CigarElement> cigar = TextCigarCodec.decode(remoteAlignment.getCigar()).getCigarElements();
				localHomologyBaseCount = Math.max(0, localBsSeq.length() - SAMRecordUtil.getStartSoftClipLength(cigar));
				if (SAMRecordUtil.getEndSoftClipLength(cigar) > 0) {
					// anchor is not aligned - something went wrong
					localHomologyBaseCount = 0;
				}
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
		return new BreakendSummary(bs.referenceIndex, bs.direction, bs.nominal + offset, bs.start + offset, bs.end + offset);
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
		int startPadding = Math.max(0, 1 - start);
		int endPadding = Math.max(0, end - refseq.getSequenceLength());
		start = Math.max(1, start);
		end = Math.min(refseq.getSequenceLength(), end);
		if (start > end) {
			// anchor is outside of contig bounds
			return StringUtils.repeat('N', end - start + 1);
		}
		byte[] bseq = lookup.getSubsequenceAt(refseq.getSequenceName(), start, end).getBases();
		if (startPadding > 0 || endPadding > 0) {
			byte[] arr = new byte[startPadding + bseq.length + endPadding];
			Arrays.fill(arr, (byte)'N');
			System.arraycopy(bseq, 0, arr, startPadding, bseq.length);
		}
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
	public static VariantContextDirectedBreakpoint annotate(ProcessingContext context, VariantContextDirectedBreakpoint bp) {
		if (!bp.isBreakendExact()) return bp;
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context, bp);
		BreakpointHomology bh = BreakpointHomology.calculate(
				context.getReference(),
				bp.getBreakendSummary().getNominalPosition(),
				bp.getUntemplatedSequence(),
				context.getVariantCallingParameters().maxBreakendHomologyLength,
				context.getVariantCallingParameters().breakendHomologyAlignmentMargin);
		int[] bounds;
		if (bp.getBreakendSummary().direction == BreakendDirection.Forward) {
			bounds = new int[] { -bh.getLocalHomologyLength(), bh.getRemoteHomologyLength() };
		} else {
			bounds = new int[] { -bh.getRemoteHomologyLength(), bh.getLocalHomologyLength() };
		}
		builder.attribute(VcfInfoAttributes.INEXACT_HOMPOS, bounds);
		return (VariantContextDirectedBreakpoint)builder.make();
	}
}
