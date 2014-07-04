package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;
import java.util.Set;

import au.edu.wehi.idsv.sam.AnomolousReadAssembly;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;

public class AssemblyBuilder {
	public AssemblyBuilder(ProcessingContext processContext) {
		this.processContext = processContext;
	}
	/**
	 * Assembly is anchored to the given genomic position 
	 * @param anchor
	 */
	public AssemblyBuilder referenceAnchor(int referenceIndex, int position) {
		this.anchorReferenceIndex = referenceIndex;
		this.anchor = position;
		return this;
	}
	public AssemblyBuilder setMateAnchor(int referenceIndex, int position) {
		// TODO: improve this by looking at all contributing mate anchor positions and taking a 95% CI
		// around the expected window size.
		// Eg: when we assemble an entire fragment size worth of mate pairs, the expected breakend window
		// size approaches zero.
		this.mateAnchorReferenceIndex = referenceIndex;
		this.mateAnchor = position;
		return this;
	}
	public AssemblyBuilder assemblyBases(byte[] bases) {
		this.bases = bases;
		return this;
	}
	public AssemblyBuilder assemblyBaseQuality(byte[] quals) {
		this.quals = quals;
		return this;
	}
	public AssemblyBuilder assembledReadCount(int readCount) {
		this.readCount = readCount;
		return this;
	}
	public AssemblyBuilder assembledBaseCount(int baseCount) {
		this.baseCount = baseCount;
		return this;
	}
	/**
	 * Length of longest soft clip contributing to this read assembly
	 * @param maxsclen
	 * @return
	 */
	public AssemblyBuilder maximumSoftClipLength(int maxsclen) {
		this.maxsclen = maxsclen;
		return this;
	}
	public AssemblyBuilder setLongestSupportingRead(int longestSupportingRead) {
		this.longestSupportingRead = longestSupportingRead;
		return this;
	}
	/**
	 * Number of assembly bases anchored to the reference.
	 * @param anchorLen
	 * @return
	 */
	public AssemblyBuilder anchorLength(int anchorLen) {
		this.anchorLen = anchorLen;
		return this;
	}
	/**
	 * Direction of assembled variant
	 */
	public AssemblyBuilder direction(BreakendDirection dir) {
		this.dir = dir;
		return this;
	}
	public AssemblyBuilder assemblerName(String name) {
		this.assembler = name;
		return this;
	}
	private ProcessingContext processContext;
	private AnomolousReadAssembly ara;
	private int anchorReferenceIndex = -1;
	private int anchor;
	private int mateAnchorReferenceIndex = -1;
	private int mateAnchor; 
	private byte[] bases;
	private byte[] quals;
	private int readCount;
	private int baseCount;
	private int maxsclen;
	private int longestAssembledRead;
	private int anchorLen;	
	private BreakendDirection dir;
	private String assembler = ""; 
	/**
	 * Generates a VCF variant from the assembly
	 * @return
	 */
	public VariantContextDirectedBreakpoint makeVariant() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(processContext);
		AnomolousReadAssembly ara = makeSAMRecord();
		
		EvidenceMetrics evidence = new EvidenceMetrics();
		evidence.set(VcfAttributes.ASSEMBLY_READS, readCount);
		evidence.set(VcfAttributes.ASSEMBLY_BASES, baseCount);
		evidence.set(VcfAttributes.ASSEMBLY_LENGTH, bases.length);
		evidence.set(VcfAttributes.ASSEMBLY_MAX_READ_SOFT_CLIP, maxsclen);
		
		builder
			.evidence(evidence)
			.attribute(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), assembler)
			.attribute(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), new String(bases, StandardCharsets.US_ASCII))
			.attribute(VcfAttributes.ASSEMBLY_QUALITY.attribute(), ara.getAssemblyBreakpointQuality())
			.attribute(VcfAttributes.ASSEMBLY_READ_LENGTH.attribute(), longestAssembledRead);
		if (anchorReferenceIndex != -1) {
			BreakendSummary summary = new BreakendSummary(anchorReferenceIndex, dir, anchor, anchor, evidence);
			builder.breakend(summary,
					dir == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBases(ara) : SAMRecordUtil.getStartSoftClipBases(ara),
					dir == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBaseQualities(ara) : SAMRecordUtil.getStartSoftClipBaseQualities(ara));
			builder.id(String.format("%s-%s:%d-%s",
					assembler,
					processContext.getDictionary().getSequence(anchorReferenceIndex).getSequenceName(),
					anchor,
					dir == BreakendDirection.Forward ? "f" : "b"));
			
			// apply filter:
			if (ara.getBreakpointLength() <= maxsclen) {
				// Just assembled one of the soft clips - not every exciting
				builder.filter("SC1"); // assembly length = 1 SC read
			}
		} else if (mateAnchorReferenceIndex != -1) {
			
			builder.id(String.format("%s-%s:%d-%s-%d",
					assembler,
					processContext.getDictionary().getSequence(anchorReferenceIndex).getSequenceName(),
					anchor,
					dir == BreakendDirection.Forward ? "f" : "b",
					new String(bases, StandardCharsets.US_ASCII).hashCode() % 10000 // try to make id unique
					));
			
			// apply filter:
			if (ara.getBreakpointLength() <= longestAssembledRead) {
				// just assembled a single read - not very exciting
				builder.filter("LN1"); // ssembly length = 1 read
			}
		} else {
			throw new IllegalStateException("Assembly not anchored to reference either directly, or through mate reads");
		}
		// apply common filter:
		if (readCount == 1) {
			builder.filter("RC1"); // read count = 1
		}
		return builder.make();
	}
	public AnomolousReadAssembly makeSAMRecord() {
		return new AnomolousReadAssembly(assembler, bases, quals, anchorLen, dir, readCount, baseCount);
	}
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
			double breakpointQuality,
			int maxsclen
			) {
		EvidenceMetrics evidence = evidence(fullAssembly, readCount, readBaseCount, maxsclen);
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
	public VariantContextAssemblyBuilder assembly(AnomolousReadAssembly ara) {
		this.ara = ara;
		return this;
	}
	private static EvidenceMetrics evidence(byte[] fullAssembly, int readCount, int readBaseCount, int maxsclen) {
		EvidenceMetrics evidence = new EvidenceMetrics();
		evidence.set(VcfAttributes.ASSEMBLY_READS, readCount);
		evidence.set(VcfAttributes.ASSEMBLY_BASES, readBaseCount);
		evidence.set(VcfAttributes.ASSEMBLY_LENGTH, fullAssembly.length);
		evidence.set(VcfAttributes.ASSEMBLY_MAX_READ_SOFT_CLIP, maxsclen);
		return evidence;
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
		EvidenceMetrics evidence = evidence(assembly, readCount, readBaseCount, 0);
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
	public VariantContextAssemblyBuilder anchor(int referenceIndex, int position, BreakendDirection direction) {
		id(String.format("%s-%s:%d-%s",
				assemblerName,
				processContext.getDictionary().getSequence(referenceIndex).getSequenceName(),
				position,
				direction == BreakendDirection.Forward ? "f" : "b"));
		return this;
	}
	*/
}
