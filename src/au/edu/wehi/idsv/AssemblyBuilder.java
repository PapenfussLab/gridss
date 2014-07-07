package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;

import au.edu.wehi.idsv.sam.AnomolousReadAssembly;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfConstants;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

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
	public AssemblyBuilder mateAnchor(int referenceIndex, int position) {
		// TODO: improve this by looking at all contributing mate anchor positions and taking a 95% CI
		// around the expected window size.
		// Eg: when we assemble an entire fragment size worth of mate pairs, the expected breakend window
		// size approaches zero.
		this.mateAnchorReferenceIndex = referenceIndex;
		this.mateAnchor = position;
		return this;
	}
	/**
	 * Adds the given mate anchor information to the assembly
	 * @param referenceIndex reference index of mate anchor
	 * @param position genomic position of mate anchor
	 * @param minDistance 95% CI for minimum number of bases between this position and the mate anchor 
	 * @param minDistance 95% CI for maximum number of bases between this position and the mate anchor
	 * @return
	 */
	public AssemblyBuilder mateAnchor(int referenceIndex, int position, int contigOffset, int minDistance, int maxDistance) {
		throw new RuntimeException("NYI");
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
	public AssemblyBuilder longestSupportingRead(int longestSupportingRead) {
		this.longestAssembledRead = longestSupportingRead;
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
	private BreakendDirection dir = null;
	private String assembler = "";
	/**
	 * Generates a VCF variant from the assembly
	 * @return
	 */
	public VariantContextDirectedBreakpoint makeVariant() {
		if (dir == null) throw new IllegalStateException("direction is required");
		if (bases == null) throw new IllegalStateException("assembly base sequence is required");
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(processContext);
		AnomolousReadAssembly ara = makeSAMRecord();
		
		EvidenceMetrics evidence = new EvidenceMetrics();
		evidence.set(VcfAttributes.ASSEMBLY_READS, readCount);
		evidence.set(VcfAttributes.ASSEMBLY_BASES, baseCount);
		evidence.set(VcfAttributes.ASSEMBLY_LENGTH, bases.length);
		builder
			.evidence(evidence)
			.attribute(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), assembler)
			.attribute(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), new String(bases, StandardCharsets.US_ASCII))
			.attribute(VcfAttributes.ASSEMBLY_READ_LENGTH.attribute(), longestAssembledRead)
			.attribute(VcfAttributes.ASSEMBLY_QUALITY.attribute(), ara.getAverageBreakpointQuality())
			.attribute(VcfAttributes.ASSEMBLY_MAX_READ_SOFT_CLIP.attribute(), maxsclen);
		if (anchorReferenceIndex != -1) {
			BreakendSummary summary = new BreakendSummary(anchorReferenceIndex, dir, anchor, anchor, evidence);
			builder.breakend(summary, ara.getBreakpointBases(), ara.getBreakpointQualities());
			builder.id(String.format("%s-%s:%d-%s",
					assembler,
					processContext.getDictionary().getSequence(anchorReferenceIndex).getSequenceName(),
					anchor,
					dir == BreakendDirection.Forward ? "f" : "b"));
			
			// apply filter:
			if (ara.getBreakpointLength() <= maxsclen) {
				// Just assembled one of the soft clips - not every exciting
				//builder.filter("SC1"); // assembly length = 1 SC read
			}
			if (ara.getBreakpointLength() == 0) {
				builder.filter(VcfFilter.ASSEMBLY_REF.filter());
			}
		} else if (mateAnchorReferenceIndex != -1) {
			// TODO: improve this by looking at all contributing mate anchor positions and taking a 95% CI
			// around the expected window size.
			// Eg: when we assemble an entire fragment size worth of mate pairs, the expected breakend window
			// size approaches zero.
			int breakpointWindow = processContext.getMetrics().getMaxFragmentSize();
			BreakendSummary summary = new BreakendSummary(
					mateAnchorReferenceIndex,
					dir,
					Math.max(mateAnchor - (dir == BreakendDirection.Backward ? breakpointWindow : 0), 1),
					Math.min(mateAnchor + (dir == BreakendDirection.Forward ? breakpointWindow : 0), processContext.getDictionary().getSequence(mateAnchorReferenceIndex).getSequenceLength()),
					evidence);
			builder.breakend(summary, ""); // no untemplated sequence as we don't have an assembly of the breakend itself
			builder.id(String.format("%s-%s:%d-%s-%d",
					assembler,
					processContext.getDictionary().getSequence(mateAnchorReferenceIndex).getSequenceName(),
					mateAnchor,
					dir == BreakendDirection.Forward ? "f" : "b",
					new String(bases, StandardCharsets.US_ASCII).hashCode() % 10000 // try to make id unique
					));
			builder.attribute(VcfSvConstants.IMPRECISE_KEY, true);
			// apply filter:
			if (ara.getBreakpointLength() <= longestAssembledRead) {
				// just assembled a single read - not very exciting
				builder.filter(VcfFilter.ASSEMBLY_TOO_SHORT.filter()); // assembly length = 1 read
			}
		} else {
			throw new IllegalStateException("Assembly not anchored to reference either directly, or through mate reads");
		}
		// apply common filter:
		if (readCount == 1) {
			builder.filter(VcfFilter.ASSEMBLY_SINGLE_READ.filter()); // read count = 1
		}
		return builder.make();
	}
	public AnomolousReadAssembly makeSAMRecord() {
		return new AnomolousReadAssembly(assembler, bases, quals, anchorLen, dir, readCount, baseCount);
	}
}
