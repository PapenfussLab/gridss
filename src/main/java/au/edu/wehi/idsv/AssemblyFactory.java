package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.Lists;

public final class AssemblyFactory {
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
	 * @param normalBaseCount number of assembly bases contributed by normal evidence sources
	 * @param tumourBaseCount number of assembly bases contributed by tumour evidence sources
	 * @return assembly evidence for the given assembly
	 */
	public static SAMRecordAssemblyEvidence createAnchored(
			ProcessingContext processContext,
			AssemblyEvidenceSource source, BreakendDirection direction,
			Set<DirectedEvidence> evidence,
			int anchorReferenceIndex, int anchorBreakendPosition, int anchoredBaseCount,
			byte[] baseCalls, byte[] baseQuals,
			int normalBaseCount, int tumourBaseCount) {
		BreakendSummary breakend = new BreakendSummary(anchorReferenceIndex, direction, anchorBreakendPosition, anchorBreakendPosition);
		return new SAMRecordAssemblyEvidence(evidence, processContext.getBasicSamHeader(), breakend, source, anchoredBaseCount, baseCalls, baseQuals, normalBaseCount, tumourBaseCount);
	}
	public static AssemblyEvidence createAnchored(
			ProcessingContext processContext, AssemblyEvidenceSource source,
			Set<DirectedEvidence> breakendSupport,
			int startAnchorReferenceIndex, int startAnchorPosition, int breakendStartOffset,
			int endAnchorReferenceIndex, int endAnchorPosition, int breakendEndOffset,
			byte[] baseCalls, byte[] baseQuals, int normalBaseCount, int tumourBaseCount) {
		throw new RuntimeException("NYI");
	}
	/**
	 * Creates an assembly whose breakpoint cannot be exactly anchored to the reference  
	 * @param processContext context
	 * @param source assembly source
	 * @param direction direction of breakend
	 * @param evidence evidence supporting the assembly breakend
	 * @param baseCalls assembly base sequence as per a positive strand read into a putative anchor
	 * @param baseQuals assembly base qualities
	 * @param normalBaseCount number of assembly bases contributed by normal evidence sources
	 * @param tumourBaseCount number of assembly bases contributed by tumour evidence sources
	 * @return assembly evidence for the given assembly
	 */
	public static SAMRecordAssemblyEvidence createUnanchored(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			Set<DirectedEvidence> evidence,
			byte[] baseCalls, byte[] baseQuals,
			int normalBaseCount, int tumourBaseCount) {
		BreakendSummary breakend = Models.calculateBreakend(processContext.getLinear(), Lists.newArrayList(evidence));
		return new SAMRecordAssemblyEvidence(evidence, processContext.getBasicSamHeader(), breakend, source, 0, baseCalls, baseQuals, normalBaseCount, tumourBaseCount);
	}
	/**
	 * Updates the given assembly to incorporate the given realignment of the assembly breakend
	 * @param processContext
	 * @return
	 */
	public static SAMRecordAssemblyEvidence incorporateRealignment(ProcessingContext processContext, SAMRecordAssemblyEvidence assembly, SAMRecord realignment) {
		if (realignment == null) return assembly;
		SAMRecordAssemblyEvidence a = (SAMRecordAssemblyEvidence)assembly;
		if (realignment.getReadUnmappedFlag() || !processContext.getRealignmentParameters().realignmentPositionUnique(realignment)) {
			// Breakend did not align well enough for us to call a breakpoint
			return new SAMRecordAssemblyEvidence(a.getEvidenceSource(), assembly, realignment);
		} else {
			return new RealignedSAMRecordAssemblyEvidence(a.getEvidenceSource(), assembly, realignment);
		}
	}
	/**
	 * Rehydrates an assembly that has been persisted as a SAMRecord
	 * @param processContext
	 * @param record assembly SAMRecord
	 * @return Assembly evidence
	 */
	public static SAMRecordAssemblyEvidence hydrate(AssemblyEvidenceSource source, SAMRecord record) {
		BreakendDirection dir = SAMRecordAssemblyEvidence.getBreakendDirection(record);
		int breakendLength = SAMRecordUtil.getSoftClipLength(record, dir);
		if (breakendLength == 0 && containsIndel(record.getCigar())) {
			return new SmallIndelSAMRecordAssemblyEvidence(source,  record);
		} else {
			return new SAMRecordAssemblyEvidence(source, record, null);
		}
	}
	private static boolean containsIndel(Cigar cigar) {
		List<CigarElement> list = cigar.getCigarElements();
		for (int i = 0; i < list.size(); i++) {
			CigarElement e = list.get(i);
			if (e.getOperator() == CigarOperator.DELETION || e.getOperator() == CigarOperator.INSERTION) {
				return true;
			}
		}
		return false;
	}
}
