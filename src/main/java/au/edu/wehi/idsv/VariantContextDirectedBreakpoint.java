package au.edu.wehi.idsv;

import org.apache.commons.lang3.NotImplementedException;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextDirectedBreakpoint extends VariantContextDirectedEvidence implements DirectedBreakpoint {
	/** */
	private static final long serialVersionUID = 1L;
	public VariantContextDirectedBreakpoint(GenomicProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(processContext, source, context);
		assert(super.getBreakendSummary() instanceof BreakpointSummary);
	}
	/**
	 * Size of event assuming the simplest explanation.
	 * @return size of breakpoint event
	 */
	public Integer getEventSize() {
		BreakendSummary low = getBreakendSummary().lowBreakend();
		BreakendSummary high = getBreakendSummary().highBreakend();
		if (low.referenceIndex != high.referenceIndex) {
			return null;
		}
		int sizeAdjustment = -1; // assume indel
		if (low.direction == high.direction) {
			// inversion
			sizeAdjustment = 0;
			// requires a second breakpoint to explain - we're not going to assume this breakpoint exists
			// as the event could be a large fold-back inversion
			return null; 
		} else if (low.direction == BreakendDirection.Backward && high.direction == BreakendDirection.Forward) {
			// tandem dup
			sizeAdjustment = 1;
		}
		return high.nominal - low.nominal + getUntemplatedSequence().length() + sizeAdjustment; 
	}
	/**
	 * Size of indel event. 
	 * @return size of indel inferred by breakpoint, null otherwise
	 */
	public Integer getIndelSize() {
		BreakendSummary low = getBreakendSummary().lowBreakend();
		BreakendSummary high = getBreakendSummary().highBreakend();
		if (low.direction == high.direction) {
			// inversion
			return null;
		} else if (low.direction == BreakendDirection.Backward && high.direction == BreakendDirection.Forward) {
			// tandem dup
			return null;
		}
		return getEventSize(); 
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return (BreakpointSummary)super.getBreakendSummary();
	}
	@Override
	public int getRemoteMapq() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public String getUntemplatedSequence() {
		return getBreakpointSequenceString();
	}
	@Override
	public String getHomologySequence() {
		return getAttributeAsString(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY, "");
	}
	@Override
	public int getHomologyAnchoredBaseCount() {
		throw new RuntimeException("NYI: need to recalculate microhomology");
	}
	public static Ordering<VariantContextDirectedBreakpoint> ByRemoteBreakendLocationStart = new Ordering<VariantContextDirectedBreakpoint>() {
		public int compare(VariantContextDirectedBreakpoint o1, VariantContextDirectedBreakpoint o2) {
			BreakpointSummary b1 = o1.getBreakendSummary();
			BreakpointSummary b2 = o2.getBreakendSummary();
			return ComparisonChain.start()
			        .compare(b1.referenceIndex2, b2.referenceIndex2)
			        .compare(b1.start2, b2.start2)
			        .compare(b1.end2, b2.end2)
			        .compare(b1.referenceIndex, b2.referenceIndex)
			        .compare(b1.start, b2.start)
			        .compare(b1.end, b2.end)
			        .result();
		  }
	};
	public static Ordering<VariantContext> ByRemoteBreakendLocationStartRaw(final GenomicProcessingContext processContext) {
		return new Ordering<VariantContext>() {
			public int compare(VariantContext o1, VariantContext o2) {
				// TODO: is this performance acceptable? This is quite an expensive compare operation
				VcfBreakendSummary b1 = new VcfBreakendSummary(processContext, o1);
				VcfBreakendSummary b2 = new VcfBreakendSummary(processContext, o2);
				int ref1 = -1;
				int ref2 = -1;
				int start1 = 0;
				int start2 = 0;
				if (b1.location instanceof BreakpointSummary) {
					ref1 = ((BreakpointSummary)b1.location).referenceIndex2;
					start1 = ((BreakpointSummary)b1.location).start2;
				}
				if (b2.location instanceof BreakpointSummary) {
					ref2 = ((BreakpointSummary)b2.location).referenceIndex2;
					start2 = ((BreakpointSummary)b2.location).start2;
				}
				int result = ComparisonChain.start()
				        .compare(ref1, ref2)
				        .compare(start1, start2)
				        .result();
				return result;
			  }
		};
	}
	@Override
	public float getBreakpointQual() {
		return (float)getPhredScaledQual();
	}
	public boolean hasBreakpointSupport(int category) {
		return getBreakpointEvidenceCountReadPair(category) > 0 ||
				getBreakpointEvidenceCountSoftClip(category) > 0 ||
				getBreakpointEvidenceCountAssemblyReadPair(category) > 0 ||
				getBreakpointEvidenceCountAssemblySoftClip(category) > 0;
	}
	public int getBreakpointReadCount() {
		return getBreakpointEvidenceCountReadPair() +
				getBreakpointEvidenceCountSoftClip() +
				getBreakpointEvidenceCountAssemblyReadPair() + 
				getBreakpointEvidenceCountAssemblySoftClip();
	}
	public int getBreakpointEvidenceCount() { return getBreakpointEvidenceCountAssembly() + getBreakpointEvidenceCountReadPair() + getBreakpointEvidenceCountSoftClip(); }
	public int getBreakpointEvidenceCountAssembly() { return getBreakpointEvidenceCountLocalAssembly() + getBreakpointEvidenceCountRemoteAssembly(); }
	public int getBreakpointEvidenceCountLocalAssembly() { return getInt(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT, 0); }
	public int getBreakpointEvidenceCountRemoteAssembly() { return getInt(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE, 0); }
	public int getBreakpointEvidenceCountCompoundAssembly() { return getInt(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_COMPOUND, 0); }
	public int getBreakpointEvidenceCountReadPair(int category) { return getInt(category, VcfFormatAttributes.BREAKPOINT_READPAIR_COUNT, 0); }
	public int getBreakpointEvidenceCountSoftClip(int category) { return getInt(category, VcfFormatAttributes.BREAKPOINT_SPLITREAD_COUNT, 0); }
	public int getBreakpointEvidenceCountAssemblyReadPair(int category) { return getInt(category, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, 0); }
	public int getBreakpointEvidenceCountAssemblySoftClip(int category) { return getInt(category, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT, 0); }
	
	public double getBreakpointEvidenceQualAssembly() { return getBreakpointEvidenceQualLocalAssembly() + getBreakpointEvidenceQualRemoteAssembly(); }
	public double getBreakpointEvidenceQualLocalAssembly() { return getDouble(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL, 0); }
	public double getBreakpointEvidenceQualRemoteAssembly() { return getDouble(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, 0); }
	public double getBreakpointEvidenceQualReadPair(int category) { return getDouble(category, VcfFormatAttributes.BREAKPOINT_READPAIR_QUAL, 0); }
	public double getBreakpointEvidenceQualSoftClip(int category) { return getBreakpointEvidenceQualLocalSoftClip(category); }
	public double getBreakpointEvidenceQualLocalSoftClip(int category) { return getDouble(category, VcfFormatAttributes.BREAKPOINT_SPLITREAD_QUAL, 0); }
	
	public int getBreakpointEvidenceCountReadPair() { return getInt(VcfInfoAttributes.BREAKPOINT_READPAIR_COUNT, 0); }
	public int getBreakpointEvidenceCountSoftClip() { return getInt(VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT, 0); }
	public int getBreakpointEvidenceCountIndel() { return getInt(VcfInfoAttributes.BREAKPOINT_INDEL_COUNT, 0); }
	public int getBreakpointEvidenceCountAssemblyReadPair() { return getInt(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, 0); }
	public int getBreakpointEvidenceCountAssemblySoftClip() { return getInt(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT, 0); }
	
	public double getBreakpointEvidenceQualReadPair() { return getDouble(VcfInfoAttributes.BREAKPOINT_READPAIR_QUAL, 0); }
	public double getBreakpointEvidenceQualSoftClip() { return getDouble(VcfInfoAttributes.BREAKPOINT_SPLITREAD_QUAL, 0); }
	public double getBreakpointEvidenceQualIndel() { return getDouble(VcfInfoAttributes.BREAKPOINT_INDEL_QUAL, 0); }
	
	public int getBreakpointSupportingFragmentCount() { return getInt(VcfInfoAttributes.BREAKPOINT_VARIANT_FRAGMENTS, 0); }
	
	@Override
	public DirectedBreakpoint asRemote() {
		throw new NotImplementedException("Not required by GRIDSS");
	}
	@Override
	public String getRemoteEvidenceID() {
		return getAttributeAsString(VcfSvConstants.PARTNER_BREAKEND_ID_KEY, null);
	}
}
