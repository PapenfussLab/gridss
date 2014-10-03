package au.edu.wehi.idsv;

import htsjdk.samtools.SAMUtils;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.util.CollectionUtil;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfAttributes.Subset;

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
		List<String> bq = getAttributeAsStringList(VcfAttributes.ASSEMBLY_BREAKEND_QUALS.attribute());
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
	protected int getIntAttrSumTN(VcfAttributes attr, Subset subset) {
		if (subset == null) subset = Subset.ALL;
		switch (subset) {
			case NORMAL:
				return getAttributeAsIntListOffset(attr.attribute(), 0, 0);
			case TUMOUR:
				return getAttributeAsIntListOffset(attr.attribute(), 1, 0);
			case ALL:
			default:
				return CollectionUtil.sumInt(getAttributeAsIntList(attr.attribute()));
		}
	}
	protected int getIntAttrMaxTN(VcfAttributes attr, Subset subset) {
		if (subset == null) subset = Subset.ALL;
		switch (subset) {
			case NORMAL:
				return getAttributeAsIntListOffset(attr.attribute(), 0, 0);
			case TUMOUR:
				return getAttributeAsIntListOffset(attr.attribute(), 1, 0);
			case ALL:
			default:
				return CollectionUtil.maxInt(getAttributeAsIntList(attr.attribute()), 0);
		}
	}
	protected double getDoubleAttrTN(VcfAttributes attr, Subset subset) {
		if (subset == null) subset = Subset.ALL;
		switch (subset) {
			case NORMAL:
				return getAttributeAsDoubleListOffset(attr.attribute(), 0, 0);
			case TUMOUR:
				return getAttributeAsDoubleListOffset(attr.attribute(), 1, 0);
			case ALL:
			default:
				return CollectionUtil.sumDouble(getAttributeAsDoubleList(attr.attribute()));
		}
	}
	public int getEvidenceCount(Subset subset) {
		return getEvidenceCountAssembly() +
				getEvidenceCountReadPair(subset) +
				getEvidenceCountSoftClip(subset);
	}
	public double getBreakendLogLikelihood(Subset subset) { return getDoubleAttrTN(VcfAttributes.LOG_LIKELIHOOD_RATIO, subset); }
	public double getBreakpointLogLikelihood(Subset subset) { return getDoubleAttrTN(VcfAttributes.LOG_LIKELIHOOD_RATIO_BREAKPOINT, subset); }
	public int getReferenceReadCount(Subset subset) { return getIntAttrSumTN(VcfAttributes.REFERENCE_COUNT_READ, subset); }
	public int getReferenceReadPairCount(Subset subset) { return getIntAttrSumTN(VcfAttributes.REFERENCE_COUNT_READPAIR, subset); }
	public int getEvidenceCountReadPair(Subset subset) { return getIntAttrSumTN(VcfAttributes.READPAIR_EVIDENCE_COUNT, subset); }
	public double getBreakendLogLikelihoodReadPair(Subset subset) { return getDoubleAttrTN(VcfAttributes.READPAIR_LOG_LIKELIHOOD_RATIO, subset); }
	public int getMappedEvidenceCountReadPair(Subset subset) { return getIntAttrSumTN(VcfAttributes.READPAIR_MAPPED_READPAIR, subset); }
	public int getMapqReadPairLocalMax(Subset subset) { return getIntAttrMaxTN(VcfAttributes.READPAIR_MAPQ_LOCAL_MAX, subset); }
	public int getMapqReadPairLocalTotal(Subset subset) { return getIntAttrSumTN(VcfAttributes.READPAIR_MAPQ_LOCAL_TOTAL, subset); }
	public int getMapqReadPairRemoteMax(Subset subset) { return getIntAttrMaxTN(VcfAttributes.READPAIR_MAPQ_REMOTE_MAX, subset); }
	public int getMapqReadPairRemoteTotal(Subset subset) { return getIntAttrSumTN(VcfAttributes.READPAIR_MAPQ_REMOTE_TOTAL, subset); }
	public int getEvidenceCountSoftClip(Subset subset) { return getIntAttrSumTN(VcfAttributes.SOFTCLIP_EVIDENCE_COUNT, subset); }
	public double getBreakendLogLikelihoodSoftClip(Subset subset) { return getDoubleAttrTN(VcfAttributes.SOFTCLIP_LOG_LIKELIHOOD_RATIO, subset); }
	public int getMappedEvidenceCountSoftClip(Subset subset) { return getIntAttrSumTN(VcfAttributes.SOFTCLIP_MAPPED, subset); }
	public int getMapqSoftClipRemoteMax(Subset subset) { return getIntAttrMaxTN(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_MAX, subset); }
	public int getMapqSoftClipRemoteTotal(Subset subset) { return getIntAttrSumTN(VcfAttributes.SOFTCLIP_MAPQ_REMOTE_TOTAL, subset); }
	public int getLengthSoftClipTotal(Subset subset) { return getIntAttrSumTN(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_TOTAL, subset); }
	public int getLengthSoftClipMax(Subset subset) { return getIntAttrMaxTN(VcfAttributes.SOFTCLIP_LENGTH_REMOTE_MAX, subset); }
	public int getEvidenceCountAssembly() { return getIntAttrSumTN(VcfAttributes.ASSEMBLY_EVIDENCE_COUNT, Subset.ALL); }
	public double getBreakendLogLikelihoodAssembly() { return getDoubleAttrTN(VcfAttributes.ASSEMBLY_LOG_LIKELIHOOD_RATIO, Subset.ALL); }
	public int getMappedEvidenceCountAssembly() { return getIntAttrSumTN(VcfAttributes.ASSEMBLY_MAPPED, Subset.ALL); }
	public int getMapqAssemblyRemoteMax() { return getIntAttrMaxTN(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_MAX, Subset.ALL); }
	public int getMapqAssemblyRemoteTotal() { return getIntAttrSumTN(VcfAttributes.ASSEMBLY_MAPQ_REMOTE_TOTAL, Subset.ALL); }
	public int getAssemblyAnchorLengthMax() { return getIntAttrMaxTN(VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX, Subset.ALL); }
	public int getAssemblyBreakendLengthMax() { return getIntAttrMaxTN(VcfAttributes.ASSEMBLY_LENGTH_REMOTE_MAX, Subset.ALL); }
	public int getAssemblyBaseCount(Subset subset) { return getIntAttrSumTN(VcfAttributes.ASSEMBLY_BASE_COUNT, subset); }
	public int getAssemblySupportCountReadPair(Subset subset) { return getIntAttrSumTN(VcfAttributes.ASSEMBLY_READPAIR_COUNT, subset); }
	public int getAssemblyReadPairLengthMax(Subset subset) { return getIntAttrMaxTN(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX, subset); }
	public int getAssemblySupportCountSoftClip(Subset subset) { return getIntAttrSumTN(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT, subset); }
	public int getAssemblySoftClipLengthTotal(Subset subset) { return getIntAttrSumTN(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, subset); }
	public int getAssemblySoftClipLengthMax(Subset subset) { return getIntAttrMaxTN(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, subset); }
	public String getAssemblerProgram() { return getAttributeAsString(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), ""); }
	public int getAssemblySupportCount(Subset subset) { return getAssemblySupportCountReadPair(subset) + getAssemblySupportCountSoftClip(subset); }
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
