package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;

import java.nio.charset.StandardCharsets;
import java.util.Iterator;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;

/**
 * VCF Breakend record
 * see Section 5.4.9 of http://samtools.github.io/hts-specs/VCFv4.2.pdf for details of breakends
 * @author Daniel Cameron
 *
 */
public class VariantContextDirectedEvidence extends IdsvVariantContext implements DirectedEvidence {
	private final VcfBreakendSummary breakend;
	//private static Log LOG = Log.getInstance(VariantContextDirectedBreakpoint.class);
	public VariantContextDirectedEvidence(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(processContext, source, context);
		this.breakend = new VcfBreakendSummary(processContext, context);
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		if (breakend.location == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return breakend.location;
	}
	@Override
	public String getEvidenceID() {
		return getID();
	}
	@Override
	public byte[] getBreakendSequence() {
		return breakend.breakpointSequence.getBytes(StandardCharsets.US_ASCII);
	}
	public String getBreakpointSequenceString() {
		if (breakend.breakpointSequence == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return breakend.breakpointSequence;
	}
	@Override
	public byte[] getBreakendQuality() {
		return null;
	}
	@Override
	public boolean isValid() {
		return breakend.location != null;
	}
	public int getReferenceReadCount(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.REFERENCE_READ_COUNT.attribute()), subset); }
	public int getReferenceReadPairCount(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.REFERENCE_READPAIR_COUNT.attribute()), subset); }
	/**
	 * Returns an iterator containing only the breakend variants from the given iterator
	 * @param context processing context
	 * @param variants input variants
	 * @return breakend variants
	 */
	public static Iterator<VariantContextDirectedEvidence> breakendIterator(final ProcessingContext context, final EvidenceSource source, final Iterator<VariantContext> variants) {
		 Iterator<IdsvVariantContext> itivc = Iterators.transform(variants, new Function<VariantContext, IdsvVariantContext>() {
			public IdsvVariantContext apply(VariantContext v) {
				return new IdsvVariantContextBuilder(context, v).source(source).make();
			}
		});
		 Iterator<VariantContextDirectedEvidence> itde = Iterators.filter(Iterators.filter(itivc, VariantContextDirectedEvidence.class), new Predicate<VariantContextDirectedEvidence>() {
			public boolean apply(VariantContextDirectedEvidence v) {
				return v.isValid();
			}
		});
		return itde;
	}
	public static Ordering<VariantContextDirectedEvidence> ByBreakendStartEnd = new Ordering<VariantContextDirectedEvidence>() {
		public int compare(VariantContextDirectedEvidence o1, VariantContextDirectedEvidence o2) {
			return BreakendSummary.ByStartEnd.compare(o1.getBreakendSummary(), o2.getBreakendSummary());
		  }
	};
	@Override
	public int getLocalMapq() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public int getLocalBaseLength() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public int getLocalMaxBaseQual() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public int getLocalTotalBaseQual() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public boolean isBreakendExact() {
		return !hasAttribute(VcfSvConstants.IMPRECISE_KEY);
	}
	@Override
	public float getBreakendQual() {
		return (float)getPhredScaledQual();
	}
	public int getBreakendEvidenceCount(EvidenceSubset subset) {
		return getBreakendEvidenceCountAssembly() +
				getBreakendEvidenceCountReadPair(subset) +
				getBreakendEvidenceCountSoftClip(subset);
	}
	public int getBreakendEvidenceCountAssembly() { return AttributeConverter.asInt(getAttribute(VcfAttributes.BREAKEND_ASSEMBLY_COUNT.attribute()), 0); }
	public int getBreakendEvidenceCountReadPair(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.BREAKEND_READPAIR_COUNT.attribute()), subset); }
	public int getBreakendEvidenceCountSoftClip(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT.attribute()), subset); }
}
