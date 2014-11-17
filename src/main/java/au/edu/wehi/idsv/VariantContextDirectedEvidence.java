package au.edu.wehi.idsv;

import htsjdk.samtools.SAMUtils;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.vcf.VcfAttributes;

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
	public String getAnchorSequenceString() {
		String ass = anchorFromAssembly();
		String alt = anchorFromAllele();
		// grab the longest anchor we can construct
		if (ass.length() > alt.length()) return ass;
		return alt;
	}
	private String anchorFromAssembly() {
		String aCons = getAssemblyConsensus();
		int anchorLength = getAssemblyAnchorLengthMax();
		if (aCons == null || aCons.length() == 0) return "";
		if (anchorLength > aCons.length()) throw new IllegalStateException(String.format("Sanity check failure: anchor length is longer than assembly consensus"));
		if (breakend.location.direction == BreakendDirection.Forward) {
			return aCons.substring(0, anchorLength);
		} else {
			return aCons.substring(aCons.length() - anchorLength);
		}
	}
	private String anchorFromAllele() {
		if (breakend.anchorSequence == null) return "";
		return breakend.anchorSequence;
	}
	@Override
	public byte[] getBreakendQuality() {
		List<String> bq = AttributeConverter.asStringList(getAttribute(VcfAttributes.ASSEMBLY_BREAKEND_QUALS.attribute()));
		if (bq.size() == 0) return null;
		try {
			return SAMUtils.fastqToPhred(URLDecoder.decode(bq.get(0), "UTF-8"));
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException(String.format("Sanity check failure: unable to decode %s", bq.get(0)), e);
		}		
	}
	@Override
	public boolean isValid() {
		return breakend.location != null;
	}
	public int getEvidenceCount(EvidenceSubset subset) {
		return getEvidenceCountAssembly() +
				getEvidenceCountReadPair(subset) +
				getEvidenceCountSoftClip(subset);
	}
	public double getBreakendLogLikelihood(EvidenceSubset subset) { return AttributeConverter.asDoubleSumTN(getAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO.attribute()), subset); }
	public double getBreakpointLogLikelihood(EvidenceSubset subset) { return AttributeConverter.asDoubleSumTN(getAttribute(VcfAttributes.LOG_LIKELIHOOD_RATIO_BREAKPOINT.attribute()), subset); }
	public int getReferenceReadCount(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.REFERENCE_COUNT_READ.attribute()), subset); }
	public int getReferenceReadPairCount(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.REFERENCE_COUNT_READPAIR.attribute()), subset); }
	public int getEvidenceCountReadPair(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.READPAIR_EVIDENCE_COUNT.attribute()), subset); }
	public double getBreakendLogLikelihoodReadPair(EvidenceSubset subset) { return AttributeConverter.asDoubleSumTN(getAttribute(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO.attribute()), subset); }
	public int getMappedEvidenceCountReadPair(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.READPAIR_MAPPED_READPAIR.attribute()), subset); }
	public int getMapqReadPairLocalMax(EvidenceSubset subset) { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.READPAIR_MAPQ_LOCAL_MAX.attribute()), subset); }
	public int getMapqReadPairLocalTotal(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.READPAIR_MAPQ_LOCAL_TOTAL.attribute()), subset); }
	public int getMapqReadPairRemoteMax(EvidenceSubset subset) { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.READPAIR_MAPQ_REMOTE_MAX.attribute()), subset); }
	public int getMapqReadPairRemoteTotal(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.READPAIR_MAPQ_REMOTE_TOTAL.attribute()), subset); }
	public int getEvidenceCountSoftClip(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.SOFTCLIP_EVIDENCE_COUNT.attribute()), subset); }
	public double getBreakendLogLikelihoodSoftClip(EvidenceSubset subset) { return AttributeConverter.asDoubleSumTN(getAttribute(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO.attribute()), subset); }
	public int getMappedEvidenceCountSoftClip(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.SOFTCLIP_MAPPED.attribute()), subset); }
	public int getMapqSoftClipRemoteMax(EvidenceSubset subset) { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_MAX.attribute()), subset); }
	public int getMapqSoftClipRemoteTotal(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_TOTAL.attribute()), subset); }
	public int getLengthSoftClipTotal(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_TOTAL.attribute()), subset); }
	public int getLengthSoftClipMax(EvidenceSubset subset) { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_MAX.attribute()), subset); }
	public int getEvidenceCountAssembly() { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.ASSEMBLY_EVIDENCE_COUNT.attribute()), EvidenceSubset.ALL); }
	public double getBreakendLogLikelihoodAssembly() { return AttributeConverter.asDoubleSumTN(getAttribute(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO.attribute()), EvidenceSubset.ALL); }
	public int getMappedEvidenceCountAssembly() { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.ASSEMBLY_MAPPED.attribute()), EvidenceSubset.ALL); }
	public int getMapqAssemblyRemoteMax() { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_MAX.attribute()), EvidenceSubset.ALL); }
	public int getMapqAssemblyRemoteTotal() { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_TOTAL.attribute()), EvidenceSubset.ALL); }
	public int getAssemblyLocalMapq() { return AttributeConverter.asInt(getAttribute(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX.attribute()), 0); }
	public int getAssemblyAnchorLengthMax() { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX.attribute()), EvidenceSubset.ALL); }
	public int getAssemblyBreakendLengthMax() { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.ASSEMBLY_LENGTH_REMOTE_MAX.attribute()), EvidenceSubset.ALL); }
	public int getAssemblyBaseCount(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.ASSEMBLY_BASE_COUNT.attribute()), subset); }
	public int getAssemblySupportCountReadPair(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.ASSEMBLY_READPAIR_COUNT.attribute()), subset); }
	public int getAssemblyReadPairLengthMax(EvidenceSubset subset) { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX.attribute()), subset); }
	public int getAssemblySupportCountSoftClip(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT.attribute()), subset); }
	public int getAssemblySoftClipLengthTotal(EvidenceSubset subset) { return AttributeConverter.asIntSumTN(getAttribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL.attribute()), subset); }
	public int getAssemblySoftClipLengthMax(EvidenceSubset subset) { return AttributeConverter.asIntMaxTN(getAttribute(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX.attribute()), subset); }
	public String getAssemblerProgram() { return getAttributeAsString(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), ""); }
	public int getAssemblySupportCount(EvidenceSubset subset) { return getAssemblySupportCountReadPair(subset) + getAssemblySupportCountSoftClip(subset); }
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
	public int getLocalBaseCount() {
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
	public boolean isSimpleAssembly() {
		return getEvidenceCountAssembly() == 1
				&& getEvidenceCountReadPair(null) == 0
				&& getEvidenceCountSoftClip(null) == 0;
	}
}
