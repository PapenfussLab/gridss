package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.util.CollectionUtil;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class AssemblyBuilder {
	private ProcessingContext processContext;
	private EvidenceSource source;
	private int anchorReferenceIndex = -1;
	private int anchor;
	private int mateAnchorReferenceIndex = -1;
	private int mateAnchor; 
	private byte[] bases;
	private byte[] quals;
	private int normalBaseCount;
	private int tumourBaseCount;
	private int anchorLen;	
	private BreakendDirection dir = null;
	private String assembler = "";
	private Set<DirectedEvidence> evidence = Sets.newHashSet();
	public AssemblyBuilder(ProcessingContext processContext, EvidenceSource source) {
		this.processContext = processContext;
		this.source = source;
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
	public AssemblyBuilder contributingEvidence(Collection<DirectedEvidence> evidence) {
		this.evidence.addAll(evidence);
		return this;
	}
	public AssemblyBuilder contributingEvidence(DirectedEvidence evidence) {
		this.evidence.add(evidence);
		return this;
	}
	public AssemblyBuilder assembledBaseCount(int normalBaseCount, int tumourBaseCount) {
		this.normalBaseCount = normalBaseCount;
		this.tumourBaseCount = tumourBaseCount;
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
	/**
	 * Generates a VCF variant from the assembly
	 * @return
	 */
	public VariantContextDirectedEvidence makeVariant() {
		if (dir == null) throw new IllegalStateException("direction is required");
		if (bases == null) throw new IllegalStateException("assembly base sequence is required");
		List<NonReferenceReadPair> rp = Lists.newArrayList(Iterables.filter(evidence, NonReferenceReadPair.class));
		List<SoftClipEvidence> sc = Lists.newArrayList(Iterables.filter(evidence, SoftClipEvidence.class));
		List<NonReferenceReadPair> rpNormal = Lists.newArrayList(Iterables.filter(rp, new Predicate<NonReferenceReadPair>() { public boolean apply(NonReferenceReadPair e) { return !((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		List<NonReferenceReadPair> rpTumour = Lists.newArrayList(Iterables.filter(rp, new Predicate<NonReferenceReadPair>() { public boolean apply(NonReferenceReadPair e) { return ((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		List<SoftClipEvidence> scNormal = Lists.newArrayList(Iterables.filter(sc, new Predicate<SoftClipEvidence>() { public boolean apply(SoftClipEvidence e) { return !((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		List<SoftClipEvidence> scTumour = Lists.newArrayList(Iterables.filter(sc, new Predicate<SoftClipEvidence>() { public boolean apply(SoftClipEvidence e) { return ((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		
		
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext);
		AnomolousReadAssembly ara = makeSAMRecord();
		
		builder.attribute(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO, Float.NaN);
		
		builder.attribute(VcfAttributes.ASSEMBLY_CONSENSUS, new String(bases, StandardCharsets.US_ASCII));
		builder.attribute(VcfAttributes.ASSEMBLY_EVIDENCE_COUNT, 1);
		builder.attribute(VcfAttributes.ASSEMBLY_BASE_COUNT, new int[] { normalBaseCount, tumourBaseCount });
		builder.attribute(VcfAttributes.ASSEMBLY_PROGRAM, assembler);
		builder.attribute(VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX, anchorLen);
		builder.attribute(VcfAttributes.ASSEMBLY_LENGTH_REMOTE_MAX, bases.length - anchorLen);
		
		builder.attribute(VcfAttributes.ASSEMBLY_READPAIR_COUNT, new int[] { rpNormal.size(), rpTumour.size() } );
		builder.attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT, new int[] { scNormal.size(), scTumour.size()} );
		
		Function<SoftClipEvidence, Integer> fscLen = new Function<SoftClipEvidence, Integer>() { public Integer apply(SoftClipEvidence e) { return e.getSoftClipLength(); } };
		Function<NonReferenceReadPair, Integer> frpReadLength = new Function<NonReferenceReadPair, Integer>() { public Integer apply(NonReferenceReadPair e) { return e.getNonReferenceRead().getReadLength(); } };
		
		int scLenN = CollectionUtil.maxInt(scNormal, fscLen, 0);
		int scLenT = CollectionUtil.maxInt(scTumour, fscLen, 0);
		int scLen = Math.max(scLenN, scLenT);
		builder.attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, new int[] { CollectionUtil.sumInt(scNormal, fscLen), CollectionUtil.sumInt(scTumour, fscLen) } );
		builder.attribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, new int[] { scLenN, scLenT } );
		int rpReadLenN = CollectionUtil.maxInt(rpNormal, frpReadLength, 0);
		int rpReadLenT = CollectionUtil.maxInt(rpTumour, frpReadLength, 0);
		int rpReadLen = Math.max(rpReadLenN, rpReadLenT);
		builder.attribute(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX, new int[] { rpReadLenN, rpReadLenT} );
		if (anchorReferenceIndex != -1) {
			BreakendSummary summary = new BreakendSummary(anchorReferenceIndex, dir, anchor, anchor);
			builder.breakend(summary, ara.getBreakpointBases(), ara.getBreakpointQualities());
			builder.id(String.format("%s-%s:%d-%s-%d",
					assembler,
					processContext.getDictionary().getSequence(anchorReferenceIndex).getSequenceName(),
					anchor,
					dir == BreakendDirection.Forward ? "f" : "b",
					Math.abs(new String(bases, StandardCharsets.US_ASCII).hashCode()) % 10000)); // approximately unique
			
			// apply filter:
			if (ara.getBreakpointLength() <= scLen) {
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
			int breakpointWindow = source.getMetrics().getMaxFragmentSize();
			BreakendSummary summary = new BreakendSummary(
					mateAnchorReferenceIndex,
					dir,
					Math.max(mateAnchor - (dir == BreakendDirection.Backward ? breakpointWindow : 0), 1),
					Math.min(mateAnchor + (dir == BreakendDirection.Forward ? breakpointWindow : 0), processContext.getDictionary().getSequence(mateAnchorReferenceIndex).getSequenceLength()));
			builder.breakend(summary, null, quals); // no untemplated sequence as we don't have an assembly of the breakend itself
			builder.id(String.format("%s-%s:%d-%s-%d",
					assembler,
					processContext.getDictionary().getSequence(mateAnchorReferenceIndex).getSequenceName(),
					mateAnchor,
					dir == BreakendDirection.Forward ? "f" : "b",
					Math.abs(new String(bases, StandardCharsets.US_ASCII).hashCode()) % 10000
					));
			builder.attribute(VcfSvConstants.IMPRECISE_KEY, true);
			// apply filter:
			if (ara.getBreakpointLength() <= rpReadLen) {
				// just assembled a single read - not very exciting
				builder.filter(VcfFilter.ASSEMBLY_TOO_SHORT.filter()); // assembly length = 1 read
			}
		} else {
			throw new IllegalStateException("Assembly not anchored to reference either directly, or through mate reads");
		}
		// apply common filter:
		if (evidence.size() == 1) {
			builder.filter(VcfFilter.ASSEMBLY_SINGLE_READ.filter()); // read count = 1
		}
		builder.source(source);
		return (VariantContextDirectedEvidence)recalculatePhredLLR(processContext, (VariantContextDirectedEvidence)builder.make());
	}
	private static VariantContextDirectedEvidence recalculatePhredLLR(ProcessingContext processContext, VariantContextDirectedEvidence assembly) {
		double llr = Models.llr(assembly);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, assembly);
		builder.attribute(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO, llr);
		builder.attribute(VcfAttributes.LOG_LIKELIHOOD_RATIO, llr);
		builder.phredScore(llr);
		return (VariantContextDirectedEvidence)builder.make();
	}
	public AnomolousReadAssembly makeSAMRecord() {
		return new AnomolousReadAssembly(assembler, bases, quals, anchorLen, dir, evidence.size(), normalBaseCount + tumourBaseCount);
	}
	/**
	 * Updates the given assembly to incorporate the given realignment of the assembly breakend
	 * @param processContext
	 * @return
	 */
	public static VariantContextDirectedEvidence incorporateRealignment(ProcessingContext processContext, VariantContextDirectedEvidence assembly, SAMRecord realignment) {
		if (realignment == null) return assembly;
		if (realignment.getReadUnmappedFlag()) {
			// No additional evidence attributes to write since realignment was assumed to have failed for breakends
			return assembly;
		}
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, assembly);
		builder.attribute(VcfAttributes.ASSEMBLY_MAPPED, 1);
		builder.attribute(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_MAX, realignment.getMappingQuality());
		builder.attribute(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_TOTAL, realignment.getMappingQuality());
		RealignedBreakpoint rbp = new RealignedBreakpoint(processContext, assembly.getBreakendSummary(), assembly.getAnchorSequenceString(), realignment);
		builder.breakpoint(rbp.getBreakpointSummary(), rbp.getInsertedSequence()); // drop base quals as we can get them from the RealignedBreakpoint if needed
		return recalculatePhredLLR(processContext, (VariantContextDirectedEvidence)builder.make());
	}
}