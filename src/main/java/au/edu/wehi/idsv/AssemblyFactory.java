package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import com.google.common.primitives.Bytes;

import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.MessageThrottler;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

public final class AssemblyFactory {
	private static final Log log = Log.getInstance(AssemblyFactory.class);
	private AssemblyFactory() { } 
	/**
	 * Creates an assembly 
	 * @param processContext context
	 * @param source assembly source
	 * @param direction direction of breakend
	 * @param evidence evidence supporting the assembly breakend
	 * @param anchorReferenceIndex contig of anchored bases 
	 * @param anchorBreakendPosition genomic position of anchored base closest breakend
	 * @param anchoredBaseCount number of anchored bases in assembly
	 * @param baseCalls assembly base sequence as per a positive strand read over the anchor
	 * @param baseQuals assembly base qualities
	 * @return assembly evidence for the given assembly
	 */
	public static SAMRecord createAnchoredBreakend(
			ProcessingContext processContext, AssemblyEvidenceSource source, AssemblyIdGenerator assemblyIdGenerator,
			BreakendDirection direction,
			List<DirectedEvidence> evidence, List<AssemblyEvidenceSupport> supportList,
			int anchorReferenceIndex, int anchorBreakendPosition, int anchoredBaseCount,
			byte[] baseCalls, byte[] baseQuals) {
		BreakendSummary breakend = new BreakendSummary(anchorReferenceIndex, direction, anchorBreakendPosition);
		assert(breakend.direction != null);
		SAMRecord r = createAssemblySAMRecord(processContext, assemblyIdGenerator, evidence, supportList, processContext.getBasicSamHeader(), source, breakend,
				breakend.direction == BreakendDirection.Forward ? anchoredBaseCount : 0,
				breakend.direction == BreakendDirection.Backward ? anchoredBaseCount : 0,
				baseCalls, baseQuals);
		return r;
	}
	public static SAMRecord createAnchoredBreakpoint(
			ProcessingContext processContext, AssemblyEvidenceSource source, AssemblyIdGenerator assemblyIdGenerator,
			List<DirectedEvidence> evidence, List<AssemblyEvidenceSupport> supportList,
			int startAnchorReferenceIndex, int startAnchorPosition, int startAnchorBaseCount,
			int endAnchorReferenceIndex, int endAnchorPosition, int endAnchorBaseCount,
			byte[] baseCalls, byte[] baseQuals) {
		BreakpointSummary bp = new BreakpointSummary(
				startAnchorReferenceIndex, BreakendDirection.Forward, startAnchorPosition,
				endAnchorReferenceIndex, BreakendDirection.Backward, endAnchorPosition);
		assert(startAnchorBaseCount > 0);
		assert(endAnchorBaseCount > 0);
		SAMRecord r = createAssemblySAMRecord(processContext, assemblyIdGenerator, evidence, supportList, processContext.getBasicSamHeader(), source, bp,
				startAnchorBaseCount,
				endAnchorBaseCount,
				baseCalls, baseQuals);
		//System.out.println(String.format("createAnchoredBreakpoint,%d,%d,%d,%d,%d,%d,%d,%d,%s",
		//		startAnchorReferenceIndex, startAnchorPosition, startAnchorBaseCount,
		//		endAnchorReferenceIndex, endAnchorPosition, endAnchorBaseCount,
		//		baseCalls.length, baseQuals.length, r == null ? null : r.getCigarString()));
		return r;
	}
	/**
	 * Creates an assembly whose breakpoint cannot be exactly anchored to the reference  
	 * @param processContext context
	 * @param source assembly source
	 * @param evidence evidence supporting the assembly breakend
	 * @param supportList
	 * @param baseCalls assembly base sequence as per a positive strand read into a putative anchor
	 * @param baseQuals assembly base qualities
	 * @return assembly evidence for the given assembly
	 */
	public static SAMRecord createUnanchoredBreakend(
			ProcessingContext processContext, AssemblyEvidenceSource source, AssemblyIdGenerator assemblyIdGenerator,
			BreakendSummary breakend,
			List<DirectedEvidence> evidence, List<AssemblyEvidenceSupport> supportList,
			byte[] baseCalls, byte[] baseQuals) {
		SAMRecord r = createAssemblySAMRecord(processContext, assemblyIdGenerator, evidence, supportList, processContext.getBasicSamHeader(), source, breakend,
				0, 0, baseCalls, baseQuals);
		return r;
	}
	private static final byte[][] PAD_BASES = new byte[][] { new byte[] {}, new byte[] { 'N' }, new byte[] { 'N', 'N' } };
	private static final byte[][] PAD_QUALS = new byte[][] { new byte[] {}, new byte[] { 0 }, new byte[] { 0, 0 } };
	private static SAMRecord createAssemblySAMRecord(
			ProcessingContext processContext,
			AssemblyIdGenerator assemblyIdGenerator,
			List<DirectedEvidence> evidence,
			List<AssemblyEvidenceSupport> supportList,
			SAMFileHeader samFileHeader, AssemblyEvidenceSource source,
			BreakendSummary breakend,
			int startAnchoredBaseCount,
			int endAnchoredBaseCount,
			byte[] baseCalls, byte[] baseQuals) {
		assert(startAnchoredBaseCount >= 0);
		assert(endAnchoredBaseCount >= 0);
		assert(startAnchoredBaseCount + endAnchoredBaseCount <= baseCalls.length);
		assert(baseCalls.length == baseQuals.length);
		assert(breakend != null);
		// null checks useful only during testing
		if (evidence == null) {
			evidence = new ArrayList<>();
		}
		if (supportList == null) {
			supportList = evidence.stream()
					.map(e -> new AssemblyEvidenceSupport(e, Range.closed(0, baseCalls.length)))
					.collect(Collectors.toList());
		}
		SAMRecord record = new SAMRecord(samFileHeader);
		// default to the minimum mapping quality that is still valid
		record.setMappingQuality((int)Math.ceil(source.getContext().getConfig().minMapq));
		record.setReferenceIndex(breakend.referenceIndex);
		record.setReadName(assemblyIdGenerator.generate(breakend, baseCalls, startAnchoredBaseCount, endAnchoredBaseCount));
		if (startAnchoredBaseCount == 0 && endAnchoredBaseCount == 0) {
			assert(!(breakend instanceof BreakpointSummary));
			// SAM spec requires at least one mapped base
			// to conform to this, we add a placeholder mismatched bases to our read
			// in the furthest anchor position
			// and represent the breakend confidence interval as an N
			// interval anchored by Xs
			record.setAlignmentStart(breakend.start);
			LinkedList<CigarElement> ce = new LinkedList<CigarElement>();
			int len = breakend.end - breakend.start + 1;
			int padBases;
			if (len <= 2) {
				ce.add(new CigarElement(len, CigarOperator.X));
				padBases = len;
			} else {
				ce.add(new CigarElement(1, CigarOperator.X));
				ce.add(new CigarElement(len - 2, CigarOperator.N));
				ce.add(new CigarElement(1, CigarOperator.X));
				padBases = 2;
			}
			if (breakend.direction == BreakendDirection.Forward) {
				ce.addLast(new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP));
				record.setCigar(new Cigar(ce));
				record.setReadBases(Bytes.concat(PAD_BASES[padBases], baseCalls));
				record.setBaseQualities(Bytes.concat(PAD_QUALS[padBases], baseQuals));
			} else {
				ce.addFirst(new CigarElement(baseCalls.length, CigarOperator.SOFT_CLIP));
				record.setCigar(new Cigar(ce));
				record.setReadBases(Bytes.concat(baseCalls, PAD_BASES[padBases]));
				record.setBaseQualities(Bytes.concat(baseQuals, PAD_QUALS[padBases]));
			}
			record.setAttribute(SamTags.UNANCHORED, 1);
		} else {
			record.setReadBases(baseCalls);
			record.setBaseQualities(baseQuals);
			if (breakend.start != breakend.end) {
				throw new IllegalArgumentException("Imprecisely anchored breakends not supported by this constructor");
			}
			if (startAnchoredBaseCount > 0 && endAnchoredBaseCount > 0) {
				// This is a breakpoint alignment spanning the entire event
				BreakpointSummary bp = (BreakpointSummary)breakend;
				record.setAlignmentStart(breakend.start - startAnchoredBaseCount + 1);
				List<CigarElement> c = new ArrayList<CigarElement>(4);
				int insSize = baseCalls.length - startAnchoredBaseCount - endAnchoredBaseCount;
				int delSize = bp.start2 - bp.start - 1;
				c.add(new CigarElement(startAnchoredBaseCount, CigarOperator.MATCH_OR_MISMATCH));
				if (insSize != 0) {
					c.add(new CigarElement(insSize, CigarOperator.INSERTION));
				}
				if (delSize != 0) {
					if (delSize < 0) {
						if (!MessageThrottler.Current.shouldSupress(log, "negative deletion")) {
							log.warn(String.format("Negative deletions not supported by SAM specs. Breakpoint assembly has been converted to breakend. "
									+ "Sanity check failure: this should not be possible for positional assembly. Breakpoint is: %s Indel is: %dI%dD", bp, insSize, delSize));
						}
						return createAssemblySAMRecord(processContext, assemblyIdGenerator, evidence, supportList, samFileHeader, source, bp.localBreakend(), startAnchoredBaseCount, 0, baseCalls, baseQuals);
					}
					c.add(new CigarElement(delSize, CigarOperator.DELETION));
				}
				c.add(new CigarElement(endAnchoredBaseCount, CigarOperator.MATCH_OR_MISMATCH));
				record.setCigar(new Cigar(c));
			} else if (startAnchoredBaseCount > 0) {
				assert(!(breakend instanceof BreakpointSummary));
				assert(breakend.direction == BreakendDirection.Forward);
				record.setAlignmentStart(breakend.start - startAnchoredBaseCount + 1);
				record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(startAnchoredBaseCount, CigarOperator.MATCH_OR_MISMATCH),
						new CigarElement(baseCalls.length - startAnchoredBaseCount, CigarOperator.SOFT_CLIP))));
			} else { // endAnchoredBaseCount > 0
				assert(!(breakend instanceof BreakpointSummary));
				assert(breakend.direction == BreakendDirection.Backward);
				record.setAlignmentStart(breakend.start);
				record.setCigar(new Cigar(ImmutableList.of(
						new CigarElement(baseCalls.length - endAnchoredBaseCount, CigarOperator.SOFT_CLIP),
						new CigarElement(endAnchoredBaseCount, CigarOperator.MATCH_OR_MISMATCH))));
			}
		}
		if (!(breakend instanceof BreakpointSummary)) {
			record.setAttribute(SamTags.ASSEMBLY_DIRECTION, breakend.direction.toChar());
		}
		int[] truncationAmount = truncateAnchorToContigBounds(processContext, record);

		AssemblyAttributes.annotateAssembly(processContext, record, evidence, supportList);
		return record;
	}
	private static int[] truncateAnchorToContigBounds(ProcessingContext processContext, SAMRecord r) {
		int end = processContext.getDictionary().getSequences().get(r.getReferenceIndex()).getSequenceLength();
		int startTruncation = 0;
		int endTruncation = 0;
		if (r.getAlignmentStart() < 1) {
			int basesToTruncate = 1 - r.getAlignmentStart();
			ArrayList<CigarElement> cigar = new ArrayList<>(r.getCigar().getCigarElements());
			int len = cigar.get(0).getLength() - basesToTruncate;
			if (len < 0) {
				if (!MessageThrottler.Current.shouldSupress(log, "truncating assembly to contig bounds")) {
					log.warn(String.format("Attempted to truncate %d bases from start of %s with CIGAR %s", basesToTruncate, r.getReadName(), r.getCigarString()));
				}
			} else {
				cigar.set(0, new CigarElement(len, cigar.get(0).getOperator()));
				r.setAlignmentStart(1);
				r.setCigar(new Cigar(cigar));
				byte[] b = r.getReadBases();
				if (b != null && b != SAMRecord.NULL_SEQUENCE) {
					r.setReadBases(Arrays.copyOfRange(b, basesToTruncate, b.length));
				}
				byte[] q = r.getBaseQualities();
				if (q != null && q != SAMRecord.NULL_QUALS) {
					r.setBaseQualities(Arrays.copyOfRange(q, basesToTruncate, q.length));
				}
				startTruncation = basesToTruncate;
			}
		}
		if (r.getAlignmentEnd() > end) {
			int basesToTruncate = r.getAlignmentEnd() - end;
			ArrayList<CigarElement> cigar = new ArrayList<>(r.getCigar().getCigarElements());
			CigarElement ce = cigar.get(cigar.size() - 1);
			int len = ce.getLength() - basesToTruncate;
			if (len < 1) {
				if (!MessageThrottler.Current.shouldSupress(log, "truncating assembly to contig bounds")) {
					log.warn(String.format("Attempted to truncate %d bases from end of %s with CIGAR %s", basesToTruncate, r.getReadName(), r.getCigarString()));
				}
			} else {
				ce = new CigarElement(len, ce.getOperator());
				cigar.set(cigar.size() - 1, ce);
				r.setCigar(new Cigar(cigar));
				byte[] b = r.getReadBases();
				if (b != null && b != SAMRecord.NULL_SEQUENCE) {
					r.setReadBases(Arrays.copyOf(b, b.length - basesToTruncate));
				}
				byte[] q = r.getBaseQualities();
				if (q != null && q != SAMRecord.NULL_QUALS) {
					r.setBaseQualities(Arrays.copyOf(q, q.length - basesToTruncate));
				}
				endTruncation = basesToTruncate;
			}
		}
		int[] truncationAmount = new int[] {startTruncation, endTruncation };
		r.setAttribute(SamTags.ASSEMBLY_ANCHOR_TRUNCATION, truncationAmount);
		return truncationAmount;
	}
}

