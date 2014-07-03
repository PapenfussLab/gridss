package au.edu.wehi.idsv;

import au.edu.wehi.idsv.vcf.VcfAttributes;

public class ReadEvidenceAssemblerUtil {
	private ReadEvidenceAssemblerUtil() { }
	public static VariantContextDirectedBreakpointBuilder breakendBuilder(
			ProcessingContext processContext,
			String assemblerName,
			int referenceIndex,
			int position,
			BreakendDirection direction,
			byte[] breakpointSequence,
			byte[] breakpointBaseQuality,
			byte[] fullAssembly,
			byte[] fullAssemblyBaseQuality,
			int readCount,
			int readBaseCount,
			double breakpointQuality
			) {
		EvidenceMetrics evidence = evidence(fullAssembly, readCount, readBaseCount);
		BreakendSummary summary = new BreakendSummary(referenceIndex, direction, position, position, evidence);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(processContext)
			.breakend(summary, breakpointSequence, breakpointBaseQuality)
			.evidence(evidence)
			.assembly(assemblerName, fullAssembly, fullAssemblyBaseQuality, breakpointQuality);
		builder.id(String.format("%s-%s:%d-%s",
				assemblerName,
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName(),
				position,
				direction == BreakendDirection.Forward ? "f" : "b"));
		return builder;
	}
	private static EvidenceMetrics evidence(byte[] fullAssembly, int readCount, int readBaseCount) {
		EvidenceMetrics evidence = new EvidenceMetrics();
		evidence.set(VcfAttributes.ASSEMBLY_READS, readCount);
		evidence.set(VcfAttributes.ASSEMBLY_BASES, readBaseCount);
		evidence.set(VcfAttributes.ASSEMBLY_LENGTH, fullAssembly.length);
		return evidence;
	}
	public static VariantContextDirectedBreakpointBuilder breakendBuilder(
			ProcessingContext processContext,
			String assemblerName,
			int referenceIndex,
			int position,
			BreakendDirection direction,
			byte[] breakpointSequence,
			byte[] fullAssembly,
			Integer readCount,
			int readBaseCount,
			double breakpointQuality
			) {
		return breakendBuilder(processContext, assemblerName, referenceIndex, position, direction, breakpointSequence, null, fullAssembly, null, readCount, readBaseCount, breakpointQuality);
	}
	public static VariantContextDirectedBreakpointBuilder mateAnchoredBuilder(
			ProcessingContext processContext,
			String assemblerName,
			int referenceIndex,
			int mateAnchorPosition,
			BreakendDirection direction,
			byte[] assembly,
			byte[] quals,
			int readCount,
			int readBaseCount,
			double assemblyQuality,
			long uniqueId
			) {
		EvidenceMetrics evidence = evidence(assembly, readCount, readBaseCount);
		// TODO: improve this by looking at all contributing mate anchor positions and taking a 95% CI
		// around the expected window size.
		// Eg: when we assemble an entire fragment size worth of mate pairs, the expected breakend window
		// size approaches zero.
		int breakpointWindow = processContext.getMetrics().getMaxFragmentSize();
		BreakendSummary summary = new BreakendSummary(
				referenceIndex,
				direction,
				mateAnchorPosition - (direction == BreakendDirection.Backward ? breakpointWindow : 0),
				mateAnchorPosition + (direction == BreakendDirection.Forward ? breakpointWindow : 0),
				evidence);
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(processContext)
			.breakend(summary, "") // no untemplated sequence since our breakend is inexact
			.evidence(evidence)
			.assembly(assemblerName, assembly, quals, assemblyQuality);
		builder.id(String.format("%s-%s:%d-%s-%d",
				assemblerName,
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName(),
				mateAnchorPosition,
				direction == BreakendDirection.Forward ? "f" : "b",
				uniqueId));
		return builder;
	}
}
