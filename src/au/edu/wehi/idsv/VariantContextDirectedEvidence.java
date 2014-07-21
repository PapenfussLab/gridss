package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfConstants;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
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
	private final byte[] breakendBaseQual;
	//private static Log LOG = Log.getInstance(VariantContextDirectedBreakpoint.class);
	public VariantContextDirectedEvidence(ProcessingContext processContext, Set<EvidenceSource> sourceSet, VariantContext context) {
		this(processContext, sourceSet, context, null);
	}
	public VariantContextDirectedEvidence(ProcessingContext processContext, Set<EvidenceSource> sourceSet, VariantContext context, byte[] breakendBaseQual) {
		super(processContext, sourceSet, context);
		this.breakendBaseQual = breakendBaseQual;
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
	public EvidenceSource getEvidenceSource() {
		return sourceSet.iterator().next();
	}
	@Override
	public byte[] getBreakendSequence() {
		return breakend.breakpointSequence.getBytes(StandardCharsets.US_ASCII);
	}
	public String getBreakpointSequenceString() {
		if (breakend.breakpointSequence == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return breakend.breakpointSequence;
	}
	public String getAnchorSequenceString() {
		if (breakend.anchorSequence == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return breakend.anchorSequence;
	}
	@Override
	public byte[] getBreakendQuality() {
		return breakendBaseQual;
	}
	@Override
	public boolean isValid() {
		return breakend.location != null;
	}
	public String getAssemblerProgram() { return getAttributeAsString(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), ""); }
	public double getAssemblyQuality() { return getAttributeAsDouble(VcfAttributes.ASSEMBLY_QUALITY.attribute(), 0); }
	public int getAssemblyMaximumSoftClipLength() { return getAttributeAsInt(VcfAttributes.ASSEMBLY_MAX_READ_SOFT_CLIP.attribute(), 0); }
	public int getAssemblyLongestSupportingRead() { return getAttributeAsInt(VcfAttributes.ASSEMBLY_READ_LENGTH.attribute(), 0); }
	public int getTotalReferenceReadDepth() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_READ_COUNT.attribute(), 0, 0); }
	public int getNormalReferenceReadDepth() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_READ_COUNT.attribute(), 1, 0); }
	public int getTumourReferenceReadDepth() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_READ_COUNT.attribute(), 2, 0); }
	public int getTotalReferenceSpanningPairCount() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_SPANNING_READ_PAIR_COUNT.attribute(), 0, 0); }
	public int getNormalReferenceSpanningPairCount() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_SPANNING_READ_PAIR_COUNT.attribute(), 1, 0); }
	public int getTumourReferenceSpanningPairCount() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_SPANNING_READ_PAIR_COUNT.attribute(), 2, 0); }
	/**
	 * Returns an iterator containing only the breakend variants from the given iterator
	 * @param context processing context
	 * @param variants input variants
	 * @return breakend variants
	 */
	public static Iterator<VariantContextDirectedEvidence> breakendIterator(final ProcessingContext context, final EvidenceSource source, final Iterator<VariantContext> variants) {
		 Iterator<VariantContextDirectedEvidence> it = Iterators.transform(variants, new Function<VariantContext, VariantContextDirectedEvidence>() {
			public VariantContextDirectedEvidence apply(VariantContext v) {
				return new IdsvVariantContextBuilder(context, v).source(source).make();
			}
		});
		it = Iterators.filter(it, new Predicate<VariantContextDirectedEvidence>() {
			public boolean apply(VariantContextDirectedEvidence v) {
				return v.isValid();
			}
		});
		return it;
	}
	public static Ordering<VariantContextDirectedEvidence> ByBreakendStartEnd = new Ordering<VariantContextDirectedEvidence>() {
		public int compare(VariantContextDirectedEvidence o1, VariantContextDirectedEvidence o2) {
			return BreakendSummary.ByStartEnd.compare(o1.getBreakendSummary(), o2.getBreakendSummary());
		  }
	};	
}
