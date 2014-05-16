package au.edu.wehi.socrates;

import au.edu.wehi.socrates.vcf.EvidenceAttributes;

public class ReadEvidenceAssemblerUtil {
	private ReadEvidenceAssemblerUtil() { }
	public static VariantContextDirectedBreakpoint create(
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
		EvidenceMetrics evidence = new EvidenceMetrics();
		evidence.set(EvidenceAttributes.ASSEMBLY_READS, readCount);
		evidence.set(EvidenceAttributes.ASSEMBLY_BASES, readBaseCount);
		evidence.set(EvidenceAttributes.ASSEMBLY_LENGTH, fullAssembly.length);
		
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
		return builder.make();
	}
	public static VariantContextDirectedBreakpoint create(
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
		return create(processContext, assemblerName, referenceIndex, position, direction, breakpointSequence, null, fullAssembly, null, readCount, readBaseCount, breakpointQuality);
	}
}
