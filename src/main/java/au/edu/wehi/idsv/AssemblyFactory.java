package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import java.util.Set;

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
		BreakendSummary breakend = Models.calculateBreakend(Lists.newArrayList(evidence));
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
}
