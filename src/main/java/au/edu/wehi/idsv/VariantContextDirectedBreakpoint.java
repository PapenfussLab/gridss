package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

import au.edu.wehi.idsv.util.CollectionUtil;
import au.edu.wehi.idsv.vcf.VcfAttributes;

import com.google.common.base.Function;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterables;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Bytes;

public class VariantContextDirectedBreakpoint extends VariantContextDirectedEvidence implements DirectedBreakpoint {
	public VariantContextDirectedBreakpoint(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(processContext, source, context);
		assert(super.getBreakendSummary() instanceof BreakpointSummary);
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
	public int getRemoteBaseLength() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public int getRemoteBaseCount() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public int getRemoteMaxBaseQual() {
		byte[] qual = getBreakendQuality();
		if (qual == null || qual.length == 0) return 0;
		List<Byte> list = Bytes.asList(qual);
		return CollectionUtil.maxInt(Iterables.transform(list, new Function<Byte, Integer>() {
			@Override
			public Integer apply(Byte arg0) {
				return (Integer)(int)(byte)arg0;
			}}), 0);
	}
	@Override
	public int getRemoteTotalBaseQual() {
		byte[] qual = getBreakendQuality();
		if (qual == null || qual.length == 0) return 0;
		int total = 0;
		for (int i = 0; i < qual.length; i++) {
			total += qual[i];
		}
		return total;
	}
	@Override
	public String getUntemplatedSequence() {
		return getBreakpointSequenceString();
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
	public static Ordering<VariantContext> ByRemoteBreakendLocationStartRaw(final ProcessingContext processContext) {
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
	public int getBreakpointEvidenceCount(EvidenceSubset subset) {
		return getBreakpointEvidenceCountAssembly() +
				getBreakpointEvidenceCountReadPair(subset) +
				getBreakpointEvidenceCountSoftClip(subset);
	}
	public int getBreakpointEvidenceCountAssembly() { return AttributeConverter.asInt(getAttribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT.attribute()), 0); }
	public int getBreakpointEvidenceCountReadPair(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.BREAKPOINT_READPAIR_COUNT.attribute()), subset); }
	public int getBreakpointEvidenceCountSoftClip(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.BREAKPOINT_SOFTCLIP_COUNT.attribute()), subset); }
	@Override
	public DirectedBreakpoint asRemote() {
		throw new RuntimeException("NYI");
	}
}
