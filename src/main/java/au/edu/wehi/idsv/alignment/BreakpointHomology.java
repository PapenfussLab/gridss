package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
import java.util.List;

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
		if (bs.direction == BreakendDirection.Backward) {
			insertedSequence = SequenceUtil.reverseComplement(insertedSequence);
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
		String localSeq = bs.getAnchorSequence(lookup, refLength);
		String localBsSeq = bs.getAnchorSequence(lookup, seqLength);
		String localRef = bs.advance(refLength).getAnchorSequence(lookup, refLength);
		String remoteSeq = SequenceUtil.reverseComplement(bs.remoteBreakend().getAnchorSequence(lookup, refLength));
		String remoteBsSeq = SequenceUtil.reverseComplement(bs.remoteBreakend().getAnchorSequence(lookup, seqLength));
		String remoteRef = SequenceUtil.reverseComplement(bs.remoteBreakend().advance(refLength).getAnchorSequence(lookup, refLength));
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
				// #344 rev-comp remote so we always have the anchor on the same side
				// This ensures that we'll choose the same alignment on both sides if there
				// are multiple equally good alignments
				SequenceUtil.reverseComplement(breakend);
				SequenceUtil.reverseComplement(remote);
				Alignment remoteAlignment = aligner.align_smith_waterman(breakend, remote);
				List<CigarElement> cigar = TextCigarCodec.decode(remoteAlignment.getCigar()).getCigarElements();
				localHomologyBaseCount = Math.max(0, localBsSeq.length() - SAMRecordUtil.getEndSoftClipLength(cigar));
				if (SAMRecordUtil.getStartSoftClipLength(cigar) > 0) {
					// anchor is not aligned - something went wrong
					localHomologyBaseCount = 0;
				}
			}
		}
		return new BreakpointHomology(localHomologyBaseCount, remoteHomologyBaseCount);
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
